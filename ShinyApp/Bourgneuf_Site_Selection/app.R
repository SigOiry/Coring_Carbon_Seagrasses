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
library(bslib)
library(shinyWidgets)

site_configs <- list(
    bourgneuf = list(
        label = "Bourgneuf Bay, France",
        shapefile = "Meadow_samplingSite.shp",
        site_column = "Sampling",
        raster_files = c(
            "tiff_meadow/bathy/Bathy_meadow_cropped.tif",
            "tiff_meadow/Freq_above_50.tif"
        ),
        bathymetry_limits = c(-3, 1),
        frequency_limits = c(0, 90),
        download_prefix = "polygons_bourgneuf",
        plot_colors = c(Mud = "darkgreen", Meadow = "darkred")
    ),
    cadiz = list(
        label = "Cadiz Bay, Spain",
        shapefile = "Meadow_samplingSite_Cadiz.shp",
        site_column = "site",
        raster_files = c(
            "tiff_meadow/bathy/bathy_cadiz_32629.tif",
            "tiff_meadow/Freq_above_50_Cadiz.tif"
        ),
        bathymetry_limits = c(-3, 1),
        frequency_limits = c(0, 100),
        download_prefix = "polygons_cadiz",
        plot_colors = c(Stable = "darkgreen", Unstable = "darkred")
    )
)

site_choices <- c(
    "Bourgneuf Bay, France" = "bourgneuf",
    "Cadiz Bay, Spain" = "cadiz"
)
raster_choices <- c("Bathymetry", "Frequence above 50%")
default_raster <- "Frequence above 50%"

load_and_stack_rasters <- function(raster_files) {
    ref <- rast(raster_files[1])
    aligned <- Map(function(path, idx) {
        current <- rast(path)
        if (idx == 1) {
            return(current)
        }
        resample(current, ref, method = "near")
    }, raster_files, seq_along(raster_files))
    stack <- rast(aligned)
    names(stack) <- tools::file_path_sans_ext(basename(raster_files))
    stack
}

load_and_prepare <- function(raster_files, target_crs = "EPSG:4326") {
    rst <- load_and_stack_rasters(raster_files)
    rst <- project(rst, target_crs)
    names(rst) <- raster_choices[seq_len(nlyr(rst))]
    rst
}

guess_local_utm_crs <- function(sf_row) {
    centroid <- sf::st_coordinates(
        sf::st_centroid(sf::st_geometry(sf::st_transform(sf_row, 4326)))
    )[1, ]
    lon <- centroid[1]
    lat <- centroid[2]
    zone <- floor((lon + 180) / 6) + 1
    epsg <- if (lat >= 0) 32600 + zone else 32700 + zone
    paste0("EPSG:", epsg)
}

rotate_sf_polygon_single <- function(sf_row, angle_deg) {
    local_crs <- guess_local_utm_crs(sf_row)
    poly_local <- sf::st_transform(sf_row, local_crs)

    rotate_matrix <- function(coords, center_xy, angle_rad) {
        shifted_x <- coords[, 1] - center_xy[1]
        shifted_y <- coords[, 2] - center_xy[2]

        rotated_x <- shifted_x * cos(angle_rad) - shifted_y * sin(angle_rad)
        rotated_y <- shifted_x * sin(angle_rad) + shifted_y * cos(angle_rad)

        coords[, 1] <- rotated_x + center_xy[1]
        coords[, 2] <- rotated_y + center_xy[2]
        coords
    }

    rotate_part <- function(part, center_xy, angle_rad) {
        if (is.matrix(part)) {
            return(rotate_matrix(part, center_xy, angle_rad))
        }
        lapply(part, rotate_part, center_xy = center_xy, angle_rad = angle_rad)
    }

    geom_type <- as.character(st_geometry_type(poly_local)[1])
    geom <- st_geometry(poly_local)[[1]]
    center_xy <- st_coordinates(st_centroid(st_geometry(poly_local)))[1, 1:2]
    angle_rad <- as.numeric(angle_deg) * pi / 180
    rotated_geom <- rotate_part(geom, center_xy, angle_rad)

    st_geometry(poly_local) <- st_sfc(
        switch(
            geom_type,
            POLYGON = st_polygon(rotated_geom),
            MULTIPOLYGON = st_multipolygon(rotated_geom),
            stop(sprintf("Unsupported geometry type for rotation: %s", geom_type))
        ),
        crs = st_crs(poly_local)
    )

    sf::st_transform(poly_local, 4326)
}

load_site_data <- function(config) {
    shp <- st_read(config$shapefile, quiet = TRUE)
    if (!config$site_column %in% names(shp)) {
        stop(sprintf("Column '%s' not found in %s", config$site_column, config$shapefile))
    }
    if (config$site_column != "Sampling") {
        shp <- shp %>% rename(Sampling = all_of(config$site_column))
    }
    shp <- shp %>%
        select(Sampling, geometry)

    rast_stack <- load_and_prepare(config$raster_files)
    shp <- st_transform(shp, crs(rast_stack))

    list(
        polygons = shp,
        rast_stack = rast_stack
    )
}

extract_metrics <- function(polygons, rast_stack) {
    sf_poly <- polygons %>% st_transform(crs(rast_stack))
    ext <- terra::extract(rast_stack, vect(sf_poly), df = TRUE)
    ext$Sampling <- sf_poly$Sampling[ext$ID]
    ext %>%
        select(-ID) %>%
        pivot_longer(cols = -Sampling, names_to = "metric", values_to = "value")
}

build_raster_domains <- function(rast_stack, bathymetry_limits) {
    domains <- lapply(names(rast_stack), function(layer_name) {
        as.numeric(terra::global(rast_stack[[layer_name]], c("min", "max"), na.rm = TRUE))
    }) |> setNames(names(rast_stack))

    domains[["Bathymetry"]] <- bathymetry_limits
    domains
}

build_plot_colors <- function(site_key, polygon_names) {
    polygon_names <- unique(as.character(polygon_names))
    config_colors <- site_configs[[site_key]]$plot_colors
    if (!all(polygon_names %in% names(config_colors))) {
        fallback <- colorFactor("Dark2", domain = polygon_names)
        missing_names <- setdiff(polygon_names, names(config_colors))
        config_colors[missing_names] <- fallback(missing_names)
    }
    config_colors[polygon_names]
}

add_polygon_edit_toolbar <- function(map) {
    draw_toolbar_args <- list(
        map = map,
        targetGroup = "Polygons",
        editOptions = editToolbarOptions(
            selectedPathOptions = selectedPathOptions(maintainColor = TRUE)
        ),
        polygonOptions = drawPolygonOptions(),
        polylineOptions = FALSE,
        rectangleOptions = FALSE,
        circleOptions = FALSE,
        markerOptions = FALSE,
        circleMarkerOptions = FALSE
    )

    if ("drag" %in% names(formals(leaflet.extras::addDrawToolbar))) {
        draw_toolbar_args$drag <- TRUE
    }

    do.call(addDrawToolbar, draw_toolbar_args)
}

selection_screen_ui <- function() {
    div(
        class = "site-selection-screen",
        div(
            class = "site-selection-card",
            h1("Interactive Sampling Site Explorer", class = "site-selection-title"),
            p(
                "Select the study site before loading the map and comparison plots.",
                class = "site-selection-text"
            ),
            pickerInput(
                "siteSelection",
                "Study site",
                choices = site_choices,
                selected = "bourgneuf",
                width = "100%",
                options = list(style = "btn-info")
            ),
            actionBttn(
                "confirmSite",
                "Open explorer",
                style = "fill",
                color = "success",
                icon = icon("map")
            )
        )
    )
}

loading_screen_ui <- function(site_key) {
    div(
        class = "site-selection-screen",
        div(
            class = "site-selection-card",
            h2(paste("Loading", site_configs[[site_key]]$label), class = "site-selection-title"),
            p("Preparing polygons, rasters, and comparison plots.", class = "site-selection-text")
        )
    )
}

main_app_ui <- function(site_key, polygon_choices) {
    fluidRow(
        class = "full-height-row",
        column(
            width = 2,
            class = "full-height-col",
            div(
                class = "sidebar-panel",
                pickerInput(
                    "studySite",
                    "Study site",
                    choices = site_choices,
                    selected = site_key,
                    width = "100%",
                    options = list(style = "btn-info")
                ),
                span(textOutput("distText"), style = "font-size:1.07em;"),
                helpText("Use the map edit tool to drag a polygon to a new position or move its vertices. To rotate, select a polygon and enter an angle below:"),
                pickerInput(
                    "selectPolygon",
                    "Polygon to rotate",
                    choices = polygon_choices,
                    selected = polygon_choices[1],
                    options = list(style = "btn-info")
                ),
                numericInputIcon("rotationAngle", "Angle (degrees)", value = 0, min = -360, max = 360, icon = icon("undo")),
                actionBttn("rotatePoly", "Rotate Polygon", style = "fill", color = "success", icon = icon("sync-alt")),
                hr(),
                prettySwitch("showImagery", "Aerial imagery", value = TRUE, status = "info", slim = TRUE),
                prettySwitch("showPolygons", "Sampling sites", value = TRUE, status = "success", slim = TRUE),
                prettySwitch("showPolyLegend", "Polygon legend", value = TRUE, status = "primary", slim = TRUE),
                tags$label("Raster layers:", style = "margin-top:15px;"),
                prettyCheckboxGroup(
                    inputId = "showRasters",
                    label = NULL,
                    choices = raster_choices,
                    selected = default_raster,
                    animation = "smooth",
                    outline = TRUE,
                    status = "primary"
                ),
                downloadBttn("downloadShp", "Download polygons (shp)", style = "fill", color = "primary")
            )
        ),
        column(
            width = 10,
            class = "full-height-col",
            div(
                id = "main-content",
                h1(
                    paste("Interactive Sampling Site Explorer -", site_configs[[site_key]]$label),
                    id = "main-title"
                ),
                fluidRow(
                    style = "height: calc(100vh - 50px);",
                    column(
                        width = 6,
                        div(style = "height:100vh;", leafletOutput("map", height = "100vh"))
                    ),
                    column(
                        width = 6,
                        div(style = "height:100vh;", plotOutput("metricsPlot", height = "100vh"))
                    )
                )
            )
        )
    )
}

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
        background: rgba(230,240,250,0.97);
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
      #main-title { font-size: 2.1em; font-weight: 700; margin-bottom: 20px; }
      .full-height-row { height: 90vh; min-height: 90vh; }
      .full-height-col { height: 90vh; min-height: 90vh; display: flex; flex-direction: column; }
      #map, #metricsPlot { height: 90vh !important; min-height: 90vh; }
      .site-selection-screen {
        min-height: 100vh;
        display: flex;
        align-items: center;
        justify-content: center;
        background: linear-gradient(145deg, #dbeaf3 0%, #eef8ef 100%);
      }
      .site-selection-card {
        width: min(92vw, 520px);
        background: rgba(255,255,255,0.96);
        border-radius: 20px;
        box-shadow: 0 10px 35px rgba(0,0,0,0.08);
        padding: 32px 28px;
      }
      .site-selection-title {
        margin-bottom: 14px;
        font-weight: 700;
      }
      .site-selection-text {
        margin-bottom: 20px;
        color: #4f6473;
      }
    "))
    ),
    uiOutput("pageContent")
)

server <- function(input, output, session) {
    metrics_df <- reactiveVal(NULL)
    rv <- reactiveValues(
        site = NULL,
        loaded_site = NULL,
        polygons = NULL,
        polygon_choices = NULL,
        rast_stack = NULL,
        palPoly = NULL,
        plot_colors = NULL,
        raster_domains = NULL,
        raster_pals = NULL,
        lastClicked = default_raster
    )

    active_config <- reactive({
        req(rv$site)
        site_configs[[rv$site]]
    })

    output$pageContent <- renderUI({
        if (is.null(rv$site)) {
            return(selection_screen_ui())
        }
        if (!identical(rv$loaded_site, rv$site) || is.null(rv$polygon_choices)) {
            return(loading_screen_ui(rv$site))
        }
        main_app_ui(rv$site, rv$polygon_choices)
    })

    observeEvent(input$confirmSite, {
        req(input$siteSelection)
        rv$loaded_site <- NULL
        rv$site <- input$siteSelection
    }, ignoreInit = TRUE)

    observeEvent(input$studySite, {
        req(input$studySite)
        if (!identical(rv$site, input$studySite)) {
            rv$loaded_site <- NULL
            rv$site <- input$studySite
        }
    }, ignoreInit = TRUE)

    observeEvent(rv$site, {
        req(rv$site)
        config <- active_config()
        site_data <- load_site_data(config)

        rv$polygons <- st_transform(site_data$polygons, 4326)
        rv$polygon_choices <- rv$polygons$Sampling
        rv$rast_stack <- site_data$rast_stack
        rv$plot_colors <- build_plot_colors(rv$site, rv$polygons$Sampling)
        rv$palPoly <- colorFactor(
            palette = unname(rv$plot_colors),
            domain = names(rv$plot_colors)
        )
        rv$raster_domains <- build_raster_domains(rv$rast_stack, config$bathymetry_limits)
        rv$raster_pals <- lapply(rv$raster_domains, function(dom) {
            colorNumeric("viridis", domain = dom, na.color = "transparent")
        }) |> setNames(names(rv$rast_stack))
        rv$lastClicked <- default_raster
        metrics_df(extract_metrics(rv$polygons, rv$rast_stack))
        rv$loaded_site <- rv$site
    }, ignoreInit = TRUE)

    observeEvent(input$showRasters, {
        if (is.null(rv$lastClicked)) {
            rv$lastClicked <- input$showRasters
            return()
        }

        newly_checked <- setdiff(input$showRasters, rv$lastClicked)
        if (length(newly_checked) == 1) {
            updatePrettyCheckboxGroup(
                session,
                inputId = "showRasters",
                selected = newly_checked
            )
        } else if (length(input$showRasters) == 0) {
            updatePrettyCheckboxGroup(
                session,
                inputId = "showRasters",
                selected = character(0)
            )
        }

        rv$lastClicked <- input$showRasters
    }, ignoreInit = TRUE)

    observeEvent(input$map_draw_edited_features, {
        req(rv$polygons, rv$rast_stack)
        geo <- input$map_draw_edited_features
        if (is.null(geo)) {
            return()
        }

        txt <- jsonlite::toJSON(geo, auto_unbox = TRUE, digits = 8)
        new_sf <- sf::st_read(txt, quiet = TRUE) %>%
            st_transform(4326)
        if (nrow(new_sf) == 0) {
            return()
        }

        modified_ids <- new_sf$layerId
        old_sf <- rv$polygons
        keep_sf <- old_sf[!old_sf$Sampling %in% modified_ids, ]

        new_sf$Sampling <- old_sf$Sampling[match(new_sf$layerId, old_sf$Sampling)]
        new_sf <- new_sf %>%
            dplyr::select(Sampling, geometry)

        rv$polygons <- bind_rows(keep_sf, new_sf) %>%
            st_as_sf() %>%
            st_transform(4326)
        metrics_df(extract_metrics(rv$polygons, rv$rast_stack))
    })

    output$map <- renderLeaflet({
        req(
            rv$polygons,
            rv$rast_stack,
            rv$palPoly,
            rv$raster_pals,
            identical(rv$loaded_site, rv$site)
        )

        polygons <- rv$polygons
        rast_stack <- rv$rast_stack
        poly_pal <- rv$palPoly
        raster_pals <- rv$raster_pals
        bbox <- st_bbox(polygons)

        map <- leaflet() %>%
            addProviderTiles("Esri.WorldImagery", group = "Imagery") %>%
            addPolygons(
                data = polygons,
                layerId = ~Sampling,
                group = "Polygons",
                color = ~poly_pal(Sampling),
                fillColor = ~poly_pal(Sampling),
                fillOpacity = 0.45,
                weight = 2,
                opacity = 1
            ) %>%
            add_polygon_edit_toolbar()

        for (layer in names(rast_stack)) {
            map <- map %>%
                addRasterImage(
                    rast_stack[[layer]],
                    colors = raster_pals[[layer]],
                    opacity = 0.7,
                    group = layer
                )
        }

        map %>%
            fitBounds(bbox[["xmin"]], bbox[["ymin"]], bbox[["xmax"]], bbox[["ymax"]])
    })

    observeEvent(input$rotatePoly, {
        req(input$selectPolygon, input$rotationAngle, rv$polygons, rv$rast_stack)
        selected <- input$selectPolygon
        idx_sel <- which(rv$polygons$Sampling == selected)
        if (length(idx_sel) != 1) {
            return()
        }

        sf_selected <- rv$polygons[idx_sel, ]
        sf_rest <- rv$polygons[-idx_sel, ]
        sf_selected_rot <- rotate_sf_polygon_single(sf_selected, as.numeric(input$rotationAngle))

        rv$polygons <- bind_rows(sf_rest, sf_selected_rot) %>%
            st_as_sf() %>%
            st_transform(4326)
        metrics_df(extract_metrics(rv$polygons, rv$rast_stack))
    })

    observe({
        req(
            rv$site,
            rv$polygons,
            rv$palPoly,
            rv$raster_domains,
            rv$raster_pals,
            identical(rv$loaded_site, rv$site)
        )
        if (is.null(input$showImagery) || is.null(input$showPolygons) ||
            is.null(input$showPolyLegend) || is.null(input$showRasters)) {
            return()
        }

        selected_rasters <- input$showRasters
        if (is.null(selected_rasters)) {
            selected_rasters <- character(0)
        }

        proxy <- leafletProxy("map")
        if (isTRUE(input$showImagery)) {
            proxy %>% showGroup("Imagery")
        } else {
            proxy %>% hideGroup("Imagery")
        }
        if (isTRUE(input$showPolygons)) {
            proxy %>% showGroup("Polygons")
        } else {
            proxy %>% hideGroup("Polygons")
        }

        for (layer in names(rv$rast_stack)) {
            if (layer %in% selected_rasters) {
                proxy %>% showGroup(layer)
            } else {
                proxy %>% hideGroup(layer)
            }
        }

        proxy %>% clearControls()
        if (isTRUE(input$showPolyLegend)) {
            proxy %>% addLegend(
                position = "bottomright",
                pal = rv$palPoly,
                values = rv$polygons$Sampling,
                title = "Sampling site"
            )
        }
        for (layer in selected_rasters) {
            proxy %>% addLegend(
                position = "topright",
                pal = rv$raster_pals[[layer]],
                values = rv$raster_domains[[layer]],
                title = layer
            )
        }
    })

    centroid_distance <- reactive({
        req(rv$polygons)
        polys <- rv$polygons
        if (nrow(polys) < 2) {
            return(NA_real_)
        }

        centroids <- st_centroid(st_geometry(polys))
        centroids <- st_transform(centroids, 3857)
        as.numeric(st_distance(centroids[1, ], centroids[2, ]))
    })

    output$distText <- renderText({
        d <- centroid_distance()
        if (is.na(d)) {
            return("Distance between boxes: N/A")
        }
        sprintf("Distance between boxes: %.1f meters", d)
    })

    output$metricsPlot <- renderPlot({
        req(metrics_df(), rv$plot_colors, rv$site)
        config <- active_config()
        df <- metrics_df()
        plot_colors <- rv$plot_colors

        p_bathy <- df %>%
            filter(metric == "Bathymetry") %>%
            ggplot(aes(x = Sampling, y = value, fill = Sampling)) +
            geom_boxplot() +
            scale_fill_manual(values = plot_colors, drop = FALSE) +
            coord_cartesian(ylim = config$bathymetry_limits) +
            labs(
                title = "Bathymetry",
                x = "Sampling site",
                y = "Bathymetry (m)"
            ) +
            theme_minimal(base_size = 18) +
            theme(
                axis.title = element_text(size = 20),
                axis.text = element_text(size = 16),
                plot.title = element_text(size = 22, face = "bold"),
                legend.position = "right"
            )

        p_freq <- df %>%
            filter(metric == "Frequence above 50%") %>%
            ggplot(aes(x = Sampling, y = value, fill = Sampling)) +
            geom_boxplot() +
            scale_fill_manual(values = plot_colors, drop = FALSE) +
            coord_cartesian(ylim = config$frequency_limits) +
            labs(
                title = "Frequence above 50%",
                x = "Sampling site",
                y = "Frequency"
            ) +
            theme_minimal(base_size = 18) +
            theme(
                axis.title = element_text(size = 20),
                axis.text = element_text(size = 16),
                plot.title = element_text(size = 22, face = "bold"),
                legend.position = "right"
            )

        p_bathy + p_freq + plot_layout(guides = "collect")
    })

    output$downloadShp <- downloadHandler(
        filename = function() {
            req(rv$site)
            paste0(site_configs[[rv$site]]$download_prefix, "_", Sys.Date(), ".zip")
        },
        content = function(file) {
            req(rv$site, rv$polygons)
            config <- site_configs[[rv$site]]
            tmpdir <- tempfile("site_polygons_")
            dir.create(tmpdir)
            shp_base <- config$download_prefix
            shp_path <- file.path(tmpdir, shp_base)

            terra::writeVector(vect(rv$polygons), paste0(shp_path, ".shp"), overwrite = TRUE)

            shp_files <- list.files(
                tmpdir,
                pattern = paste0("^", shp_base, ".*\\.(shp|shx|dbf|prj|cpg)$"),
                full.names = TRUE
            )
            zip::zipr(zipfile = file, files = shp_files, root = tmpdir)
        },
        contentType = "application/zip"
    )
}

shinyApp(ui, server)
