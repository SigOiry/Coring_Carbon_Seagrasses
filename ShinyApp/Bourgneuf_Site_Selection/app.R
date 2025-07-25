library(shiny)
library(leaflet)
library(leaflet.extras)
library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(patchwork)
library(spatialEco)

# library(shiny)
library(bslib)
library(shinyWidgets)
# library(leaflet)

#— helper to load + align rasters, then name them ——

load_and_stack_rasters <- function(raster_files) {
    ref <- rast(raster_files[2])  # Use first raster as reference
    aligned <- lapply(raster_files, function(f) {
        r <- rast(f)
        resample(r, ref, method = "near")
    })
    stack <- rast(aligned)
    names(stack) <- tools::file_path_sans_ext(basename(raster_files))
    return(stack)
}

load_and_prepare <- function(raster_files, target_crs="EPSG:4326"){
    rst <- load_and_stack_rasters(raster_files)
    rst <- project(rst, target_crs)
    nm  <- c("Bathymetry","Frequence above 50%")
    names(rst) <- nm[seq_len(nlyr(rst))]
    rst
}

rotate_sf_polygon_single <- function(sf_row, angle_deg) {
    # 1. Project to UTM zone 30N
    poly_utm <- sf::st_transform(sf_row, 32630)
    # 2. Rotate in UTM (rotation in *degrees*)
    poly_rot <- spatialEco::rotate.polygon(
        poly_utm,
        angle = as.numeric(angle_deg),
        sp = FALSE,
        anchor = "center"
    )
    # 3. Add back the CRS lost during rotation
    sf::st_crs(poly_rot) <- sf::st_crs(poly_utm)
    # 4. Project back to WGS84
    poly_rot_wgs84 <- sf::st_transform(poly_rot, 4326)
    poly_rot_wgs84
}

update_metrics <- function() {
    sf_poly <- rv$polygons %>% st_transform(crs(rast_stack))
    ext <- terra::extract(rast_stack, vect(sf_poly), df = TRUE)
    ext$Sampling <- sf_poly$Sampling[ext$ID]
    df <- ext %>%
        select(-ID) %>%
        pivot_longer(cols = -Sampling, names_to = "metric", values_to = "value")
    df
}



raster_choices <- c("Bathymetry", "Frequence above 50%")

ui <- fluidPage(
    theme = bs_theme(
        version = 5,
        bootswatch = "minty",
        base_font = font_google("Roboto"),
        heading_font = font_google("Montserrat"),
        font_scale = 1.1
    ),
    tags$head(
        tags$style(HTML("
      html, body, .container-fluid, .row, #main-content { height: 100%; min-height: 100vh; }
      body { margin: 0; padding: 0; }
      .sidebar-panel {
        background: rgba(230,240,250,0.97); /* soft blue */
        border-radius: 16px;
        box-shadow: 0 4px 24px rgba(0,0,0,0.07);
        margin-top: 0;
        margin-bottom: 0;
        margin-left: 0;
        padding: 28px 18px 28px 18px;
        min-width: 280px;
        max-width: 350px;
        height: 90vh;
      }
      .leaflet-container, #metricsPlot { border-radius: 18px; }
      #main-title { font-size: 2.1em; font-weight: 700; margin-bottom: 20px;}
      .full-height-row { height: 90vh; min-height: 90vh; }
      .full-height-col { height: 90vh; min-height: 90vh; display: flex; flex-direction: column; }
      #map, #metricsPlot { height: 90vh !important; min-height: 90vh; }
    "))
    ),
    fluidRow(
        class = "full-height-row",
        column(
            width = 2,           # <--- make this 2 or even 1 for a narrow sidebar
            class = "full-height-col",
            div(class = "sidebar-panel",
                span(textOutput("distText"), style="font-size:1.07em;"),
                helpText("Use the map controls to edit polygons. To rotate, select a polygon and enter an angle below:"),
                pickerInput("selectPolygon", "Polygon to rotate", choices = c("Mud","Meadow"), options = list(style = "btn-info")),
                numericInputIcon("rotationAngle", "Angle (degrees)", value = 0, min = -360, max = 360, icon = icon("undo")),
                actionBttn("rotatePoly", "Rotate Polygon", style = "fill", color = "success", icon = icon("sync-alt")),
                hr(),
                prettySwitch("showImagery", "Aerial imagery", value = TRUE, status = "info", slim = TRUE),
                prettySwitch("showPolygons", "Sampling sites", value = TRUE, status = "success", slim = TRUE),
                prettySwitch("showPolyLegend", "Polygon legend", value = TRUE, status = "primary", slim = TRUE),
                tags$label("Raster layers:", style = "margin-top:15px;"),
                prettyCheckboxGroup(
                    inputId = "showRasters",
                    label   = NULL,
                    choices = raster_choices,
                    selected = NULL, 
                    animation = "smooth",
                    outline = TRUE,
                    status = "primary"
                ),
                downloadBttn("downloadShp", "Download polygons (shp)", style = "fill", color = "primary")
            )
        ),
        column(
            width = 10,    # <--- fill the rest of the space!
            class = "full-height-col",
            div(
                id = "main-content",
                h1("Interactive Sampling Site Explorer", id = "main-title"),
                fluidRow(
                    style = "height: calc(100vh - 50px);",
                    column(
                        width = 6,
                        div(style = "height:100vh;",
                            leafletOutput("map", height = "100vh")
                        )
                    ),
                    column(
                        width = 6,
                        div(style = "height:100vh;",
                            plotOutput("metricsPlot", height = "100vh")
                        )
                    )
                )
            )
        )
    )
)



server <- function(input, output, session) {
    
    metrics_df <- reactiveVal(NULL)
    #— load data —
    rasters    <- list.files("tiff_meadow", "\\.tif$", full.names = TRUE, recursive = TRUE)
    rast_stack <- load_and_prepare(rasters)
    shp        <- st_read("Meadow_samplingSite.shp") %>% st_transform(crs(rast_stack))
    rv         <- reactiveValues(polygons = shp)
    
    # polygon palette
    palPoly <- colorFactor("Set1", domain = shp$Sampling)
    
    # raster domains & palettes
    raster_domains <- lapply(names(rast_stack), function(nm) {
        as.numeric(terra::global(rast_stack[[nm]], c("min","max"), na.rm=TRUE))
    }) |> setNames(names(rast_stack))
    
    raster_domains[["Bathymetry"]] <- c(-3, 1)  # Set min=-3, max=1
    
    raster_pals <- lapply(raster_domains, function(dom) {
        colorNumeric("viridis", domain = dom, na.color = "transparent")
    }) |> setNames(names(rast_stack))
    
    # populate raster checkboxes once
    observe({
        updateCheckboxGroupInput(
            session, "showRasters",
            choices  = names(rast_stack),
            selected = "Frequence above 50%"
        )
    })
    
    observeEvent(input$map_draw_edited_features, {
        geo <- input$map_draw_edited_features
        if (is.null(geo)) return()
        txt   <- jsonlite::toJSON(geo, auto_unbox = TRUE, digits = 8)
        newSF <- sf::st_read(txt, quiet = TRUE) %>% st_transform(crs(rast_stack))
        if (nrow(newSF) == 0) return()
        if (nrow(newSF) == nrow(rv$polygons)) {
            newSF <- st_set_geometry(rv$polygons, st_geometry(newSF))
            rv$polygons <- newSF
        } else {
            for (i in seq_len(nrow(newSF))) {
                dists <- st_distance(st_centroid(newSF[i,]), st_centroid(rv$polygons))
                idx <- which.min(dists)
                st_geometry(rv$polygons)[idx] <- st_geometry(newSF)[i]
            }
        }
        rv$polygons <- sf::st_transform(rv$polygons, 4326)
        sf_poly <- rv$polygons %>% st_transform(crs(rast_stack))
        # extract pixels
        ext <- terra::extract(rast_stack, vect(sf_poly), df = TRUE)
        ext$Sampling <- sf_poly$Sampling[ext$ID]
        # tidy
        df <- ext %>%
            select(-ID) %>%
            pivot_longer(cols = -Sampling, names_to = "metric", values_to = "value")
        metrics_df(df) # <- add this line!
    })
    
    # render the base map
    output$map <- renderLeaflet({
        m <- leaflet() %>%
            addProviderTiles("Esri.WorldImagery", group = "Imagery") %>%
            addPolygons(
                data    = rv$polygons,
                layerId = ~Sampling,
                group   = "Polygons",
                color   = ~palPoly(Sampling),
                fill    = ~palPoly(Sampling), 
                opacity = 1
            ) %>%
            addDrawToolbar(
                targetGroup    = "Polygons",
                editOptions    = editToolbarOptions(),
                polygonOptions = drawPolygonOptions(),
                polylineOptions     = FALSE,
                rectangleOptions    = FALSE,
                circleOptions       = FALSE,
                markerOptions       = FALSE,
                circleMarkerOptions = FALSE
            )
        
        # add each raster
        for (layer in names(rast_stack)) {
            m <- m %>% addRasterImage(
                rast_stack[[layer]],
                colors  = raster_pals[[layer]],
                opacity = 0.7,
                group   = layer
            )
        }
        
        m
    })
    
    observeEvent(input$rotatePoly, {
        req(input$selectPolygon, input$rotationAngle)
        selected <- input$selectPolygon
        idx_sel <- which(rv$polygons$Sampling == selected)
        if (length(idx_sel) != 1) return()
        sf_selected <- rv$polygons[idx_sel, ]
        sf_rest <- rv$polygons[-idx_sel, ]
        sf_selected_rot <- rotate_sf_polygon_single(sf_selected, as.numeric(input$rotationAngle))
        rv$polygons <- dplyr::bind_rows(sf_rest, sf_selected_rot)
        rv$polygons <- sf::st_transform(rv$polygons, 4326)
        sf_poly <- rv$polygons %>% st_transform(crs(rast_stack))
        # extract pixels
        ext <- terra::extract(rast_stack, vect(sf_poly), df = TRUE)
        ext$Sampling <- sf_poly$Sampling[ext$ID]
        # tidy
        df <- ext %>%
            select(-ID) %>%
            pivot_longer(cols = -Sampling, names_to = "metric", values_to = "value")
        metrics_df(df)  # <- add this line!
    })
    
    # show/hide & legends
    observe({
        proxy <- leafletProxy("map")
        if (input$showImagery)  proxy %>% showGroup("Imagery")  else proxy %>% hideGroup("Imagery")
        if (input$showPolygons) proxy %>% showGroup("Polygons") else proxy %>% hideGroup("Polygons")
        for (layer in names(rast_stack)) {
            if (layer %in% input$showRasters) proxy %>% showGroup(layer)
            else                               proxy %>% hideGroup(layer)
        }
        proxy %>% clearControls()
        if (input$showPolyLegend) {
            proxy %>% addLegend(
                position = "bottomright",
                pal      = palPoly,
                values   = rv$polygons$Sampling,
                title    = "Sampling site"
            )
        }
        for (layer in input$showRasters) {
            proxy %>% addLegend(
                position = "topright",
                pal      = raster_pals[[layer]],
                values   = raster_domains[[layer]],
                title    = layer
            )
        }
        proxy 
    })
    
    # — NEW: extract & plot inside observeEvent on the button —
    #     so every click truly recomputes from fresh geometry
    # metrics_df <- reactiveVal(NULL)
    observeEvent(input$compute_metrics, {
        # grab latest polygons & reproject
        sf_poly <- rv$polygons %>% st_transform(crs(rast_stack))
        # extract pixels
        ext <- terra::extract(rast_stack, vect(sf_poly), df = TRUE)
        ext$Sampling <- sf_poly$Sampling[ext$ID]
        # tidy
        df <- ext %>%
            select(-ID) %>%
            pivot_longer(cols = -Sampling, names_to = "metric", values_to = "value")
        metrics_df(df)
    })
    
    centroid_distance <- reactive({
        polys <- rv$polygons
        if (nrow(polys) < 2) return(NA)
        # Get centroids of each polygon
        cents <- st_centroid(polys)
        # Get distance between first two centroids
        dist <- as.numeric(st_distance(cents[1,], cents[2,]))
        # If CRS is not metric, transform to EPSG:3857 for meters
        if (!grepl("units=m", sf::st_crs(polys)$wkt, ignore.case = TRUE)) {
            cents <- st_transform(cents, 3857)
            dist <- as.numeric(st_distance(cents[1,], cents[2,]))
        }
        dist
    })
    
    output$distText <- renderText({
        d <- centroid_distance()
        if (is.na(d)) return("Distance between boxes: N/A")
        sprintf("Distance between boxes: %.1f meters", d)
    })
    
    # render the plot from the df stored by observeEvent
    output$metricsPlot <- renderPlot({
        req(metrics_df())
        df <- metrics_df()
        # Bathymetry plot (y and color scales fixed from -3 to 1)
        p_bathy <- df %>%
            filter(metric == "Bathymetry") %>%
            ggplot(aes(x = Sampling, y = value, fill = Sampling)) +
            geom_boxplot() +
            scale_fill_manual(values = c("Mud" = "darkgreen", "Meadow" = "darkred"))+
            coord_cartesian(ylim = c(-3, 1)) +
            labs(
                title = "Bathymetry",
                x = "Sampling site",
                y = "Bathymetry (m)"
            ) +
            theme_minimal(base_size = 18) +
            theme(
                axis.title      = element_text(size = 20),
                axis.text       = element_text(size = 16),
                plot.title      = element_text(size = 22, face = "bold"),
                legend.position = "right"
            )
        
        # Frequence above 50% plot (default color scale)
        p_freq <- df %>%
            filter(metric == "Frequence above 50%") %>%
            ggplot(aes(x = Sampling, y = value, fill = Sampling)) +
            geom_boxplot() +
            scale_fill_manual(values = c("Mud" = "darkgreen", "Meadow" = "darkred"))+
            labs(
                title = "Frequence above 50%",
                x = "Sampling site",
                y = "Frequency"
            ) +
            theme_minimal(base_size = 18) +
            ylim(c(0,90))+
            theme(
                axis.title      = element_text(size = 20),
                axis.text       = element_text(size = 16),
                plot.title      = element_text(size = 22, face = "bold"),
                legend.position = "right"
            )
        
        # Combine with patchwork
        combined_plot <- p_bathy + p_freq +
            plot_layout(guides = "collect") # Shared legends if same palette, or "keep" for two legends
        
        combined_plot
    })
    
    # ** trigger an initial plot on app start **
    # simulate a click so metrics_df() is populated immediately
    observeEvent(TRUE, {
        # small delay to let everything initialize
        invalidateLater(100, session)
        isolate({
            metrics_df(
                {
                    sf_poly <- rv$polygons %>% st_transform(crs(rast_stack))
                    ext <- terra::extract(rast_stack, vect(sf_poly), df = TRUE)
                    ext$Sampling <- sf_poly$Sampling[ext$ID]
                    ext %>%
                        select(-ID) %>%
                        pivot_longer(cols = -Sampling, names_to = "metric", values_to = "value")
                }
            )
        })
    }, once = TRUE)
    
    output$downloadShp <- downloadHandler(
        filename = function() {
            paste0("polygons_BB_", Sys.Date(), ".zip")
        },
        content = function(file) {
            # Create a temp directory
            tmpdir <- tempdir()
            # Path for shapefile (without extension)
            shp_path <- file.path(tmpdir, "polygons_BB")
            # Write the shapefile using sf::st_write (it writes all necessary files: .shp, .shx, .dbf, etc.)
            terra::writeVector(vect(rv$polygons), paste0(shp_path, ".shp"), overwrite = TRUE)
            # Find all files with the same name prefix
            shp_files <- list.files(tmpdir, pattern = "polygons_BB.*(shp|shx|dbf|prj|cpg)$", full.names = TRUE)
            # Zip them
            zip::zipr(zipfile = file, files = shp_files, root = tmpdir)
        },
        contentType = "application/zip"
    )
}

shinyApp(ui, server)
