library(shiny)
library(visNetwork)
library(dplyr)
library(readr)
library(bslib)
library(base64enc)

setwd("/Users/Dasha/work/Sardinia/W4H/olink/batch12/results12/intensity_shared_prots_261125/network")
edges_data <- read.delim("network_data/network.edges.txt", sep = "\t", check.names = F, as.is = T)  
nodes_data <- read.delim("network_data/network.nodes.txt", sep = "\t", check.names = F, as.is = T) 
annot <- read.delim("Su_MR.subset_cut.txt", sep = "\t", check.names = F, as.is = T)

node_counts <- table(c(edges_data$prot, edges_data$pheno))
node_sel <- names(node_counts[node_counts > 1])
nodes_data$nodes_to_select <- ifelse(nodes_data$feature %in% node_sel, TRUE, FALSE)

# Process edges data
edges_vis <- edges_data %>%
  mutate(
    from = prot,
    to = pheno,
    strong_assoc = if_else(abs(estimate) > 0.15, TRUE, FALSE),
    color = case_when(
      estimate > 0 ~ "#FF0000",
      estimate < 0 ~ "#3498DB",
    ),
    arrows = ifelse(has_direction, "to", ""),
    width = 1 ,
    smooth = F,
    title = paste0(prot, " → ", pheno, "<br>Est: ", round(estimate, 3)),
  )

# Process nodes data
nodes_vis <- nodes_data %>%
  filter(feature %in% c(edges_vis$from, edges_vis$to)) %>%
  mutate(
    id = feature,
    label = feature,  
    color = case_when(
      type == "protein" ~ "#97C2FC",
      type == "phenotype" ~ "#FFB6C1",
      type == "hormone" ~ "#98FB98"
    ),
    node_with_two_edges = nodes_to_select,
    shape = 'circle',
    size = 25,
    borderWidth = 2,
    borderWidthSelected = 4,
    title = paste0("Type: ", type, "<br>ID: ", feature)
  )

edges_with_types <- edges_vis %>%
  left_join(nodes_data %>% select(feature, type_from = type), by = c("from" = "feature")) %>%
  left_join(nodes_data %>% select(feature, type_to = type), by = c("to" = "feature")) %>%
  mutate(
    edge_type_cat = paste(pmin(type_from, type_to), "-", pmax(type_from, type_to))
  )

# MR disease annot
annot_flt <- annot %>% filter(Assay %in% nodes_vis$feature)

disease_edges_all <- data.frame(
  from = annot_flt$Assay,
  to = annot_flt$MR_outcomes,
  color = case_when(
    annot_flt$Beta > 0 ~ "#FF0000",
    annot_flt$Beta < 0 ~ "#3498DB",
  ),      
  dashes = TRUE,
  width = 1.5,
  arrows = 'to',
  smooth = FALSE,
  title = paste("Protein:", annot_flt$Assay, "<br>Disease:", annot_flt$MR_outcomes),
  edge_type_cat = "protein - disease",
  strong_assoc = FALSE
)

unique_diseases_all <- unique(annot_flt$MR_outcomes)

disease_nodes_all <- data.frame(
  id = unique_diseases_all,
  feature = unique_diseases_all,
  label = unique_diseases_all,
  title = paste("Disease (MR Outcome):", unique_diseases_all),
  type = "disease",
  color = "#FFA500",
  shape = "box",
  size = 40,
  borderWidth = 2,
  borderWidthSelected = 4,
  widthConstraint = 150,
  node_with_two_edges = NA
)

# -------------------------------------------------- 
logo_base64 <- dataURI(file = "logo_2_1.png", mime = "")

ui <- page_sidebar(
  fillable = TRUE,
  sidebar = sidebar(
    title = div(
      img(src = logo_base64, 
          style = "width: 100px; height: auto; margin-top: -30px; margin-bottom: 5px;margin-left: 25px")
    ),
    selectizeInput("multi_node_select", 
                   "Select one or mulitple nodes:", 
                   choices = NULL, 
                   multiple = TRUE,
                   options = list(placeholder = 'Type or select nodes...')),
    helpText(
      tags$span(style = "font-size: 11px;", 
                "Select nodes from the dropdown menu or by clicking one or multiple nodes on the graph. Click empty space to reset."
      )
    ),
    hr(),
    checkboxGroupInput("edge_filters", 
                       "Show edge types:",
                       choices = c("hormone - hormone",
                                   "phenotype - phenotype",
                                   "hormone - phenotype", 
                                   "hormone - protein", 
                                   "phenotype - protein"),
                       selected = c("hormone - hormone",
                                    "phenotype - phenotype",
                                    "hormone - phenotype", 
                                    "hormone - protein", 
                                    "phenotype - protein")),
    hr(),
    checkboxInput("show_diseases", "Show associated diseases (based on MR)", value = FALSE),
    checkboxInput("filter_strong", "Show only associations with |estimate| > 0.15", value = FALSE),
    checkboxInput("filter_two_edges", "Show only nodes with multiple connections", value = FALSE)
    
  ),
  card(
    card_header(
      class = "d-flex justify-content-between align-items-center",
      
      span("Associations and causal relationships among  sex hormones, plasma proteins and CMD phenotypes from Zhernakova et al., submitted."),
      
      actionButton("btn_show_legend", "Method & Legend", 
                   icon = icon("info-circle"), 
                   class = "btn-primary btn-sm") 
    ),
    card_body(
      class = "p-0", 
      visNetworkOutput("network_plot", height = "100%")
    )
  )
)

# -------------------------------------------------- 

server <- function(input, output, session) {
  observeEvent(input$btn_show_legend, {
    showModal(modalDialog(
      title = "Network Methods & Legend",
      size = "l", 
      easyClose = TRUE, 
      
      h4("Method Description"),
      p("This network visualizes significant associations between plasma levels of sex hormones, proteins, and CMD-related phenotypes."),
      p("The analsyses were performed in the longitudinal Women4Health cohort, which consists of 159 women followed during a natural menstrual cycle."),
      p("Causal relationships between proteins and hormones/phenotypes were estimated using longitudinal cross-lagged panel models; protein - disease links were taken from a large Mendelian Randomization (MR) study"),
      
      hr(),
      
      h4("Legend"),
      fluidRow(
        column(6,
               h5("Nodes (Features)"),
               tags$div(style = "display: flex; align-items: center; margin-bottom: 5px;",
                        tags$span(style = "width: 15px; height: 15px; background-color: #97C2FC; border-radius: 50%; display: inline-block; margin-right: 10px; border: 2px solid #2B7CE9;"),
                        "Protein"),
               tags$div(style = "display: flex; align-items: center; margin-bottom: 5px;",
                        tags$span(style = "width: 15px; height: 15px; background-color: #98FB98; border-radius: 50%; display: inline-block; margin-right: 10px; border: 2px solid #2B7CE9;"),
                        "Hormone"),
               tags$div(style = "display: flex; align-items: center; margin-bottom: 5px;",
                        tags$span(style = "width: 15px; height: 15px; background-color: #FFB6C1; border-radius: 50%; display: inline-block; margin-right: 10px; border: 2px solid #2B7CE9;"),
                        "Phenotype"),
               tags$div(style = "display: flex; align-items: center; margin-bottom: 5px;",
                        tags$span(style = "width: 15px; height: 15px; background-color: #FFA500; display: inline-block; margin-right: 10px; border: 2px solid #2B7CE9;"),
                        "Disease (MR Outcome)")
        ),
        column(6,
               h5("Edges (Associations)"),
               tags$div(style = "display: flex; align-items: center; margin-bottom: 5px;",
                        tags$span(style = "width: 30px; height: 3px; background-color: #FF0000; display: inline-block; margin-right: 10px;"),
                        "Positive Association"),
               tags$div(style = "display: flex; align-items: center; margin-bottom: 5px;",
                        tags$span(style = "width: 30px; height: 3px; background-color: #3498DB; display: inline-block; margin-right: 10px;"),
                        "Negative Association"),
               tags$div(style = "display: flex; align-items: center; margin-bottom: 5px;",
                        tags$span(style = "width: 30px; height: 3px; border-top: 3px dashed #666; display: inline-block; margin-right: 10px;"),
                        "Disease Association (MR)")
        )
      ),
      footer = modalButton("Close")
    ))
  })
  
  # --- 1. Reactive Data ---
  graph_data <- reactive({
    
    filtered_edges <- edges_with_types %>%
      filter(edge_type_cat %in% input$edge_filters)
    
    if (isTRUE(input$filter_strong)) {
      filtered_edges <- filtered_edges %>%
        filter(strong_assoc == TRUE)
    }
    
    active_node_ids <- unique(c(filtered_edges$from, filtered_edges$to))
    filtered_nodes <- nodes_vis %>% filter(id %in% active_node_ids)
    
    if (isTRUE(input$filter_two_edges)) {
      filtered_nodes <- filtered_nodes %>%
        filter(node_with_two_edges == TRUE)
      
      current_degrees <- table(c(filtered_edges$from, filtered_edges$to))
      nodes_with_2_plus_edges <- names(current_degrees[current_degrees >= 2])
      
      filtered_nodes <- filtered_nodes %>%
        filter(id %in% nodes_with_2_plus_edges)
      
      filtered_edges <- filtered_edges %>%
        filter(from %in% filtered_nodes$id & to %in% filtered_nodes$id)
      
      final_active_ids <- unique(c(filtered_edges$from, filtered_edges$to))
      filtered_nodes <- filtered_nodes %>% filter(id %in% final_active_ids)
    }
    
    if (isTRUE(input$show_diseases)) {
      visible_proteins <- filtered_nodes$id
      disease_edges <- disease_edges_all %>% filter(from %in% visible_proteins)
      
      if (nrow(disease_edges) > 0) {
        disease_nodes <- disease_nodes_all %>% filter(id %in% disease_edges$to)
        filtered_edges <- bind_rows(filtered_edges, disease_edges)
        filtered_nodes <- bind_rows(filtered_nodes, disease_nodes)
      }
    }
    
    filtered_edges <- filtered_edges %>% mutate(id = paste0("e", row_number()))
    filtered_nodes <- filtered_nodes %>% arrange(id)
    
    list(nodes = filtered_nodes, edges = filtered_edges)
  })
  
  # --- 2. Update Dropdown Choices ---
  observe({
    dat <- graph_data()
    updateSelectizeInput(session, "multi_node_select", 
                         choices = sort(unique(dat$nodes$id)),
                         server = TRUE)
  })
  
  # --- 3. Render Network ---
  output$network_plot <- renderVisNetwork({
    dat <- graph_data()
    
    visNetwork(nodes = dat$nodes, edges = dat$edges) %>%
      visNodes(
        shape = "circle",
        widthConstraint = list(minimum = 75, maximum = 75), 
        heightConstraint = list(minimum = 75, valign = "middle"),
        font = list(size = 16, color = "#000000", face = "arial", strokeWidth = 0),
        opacity = 1,
        scaling = list(
          label = list(enabled = TRUE, min = 20, max = 50, maxVisible = 10000)
        )
      ) %>%
      visOptions(
        highlightNearest = FALSE,
        nodesIdSelection = FALSE 
      ) %>%
      visPhysics(
        solver = "forceAtlas2Based",
        forceAtlas2Based = list(
          gravitationalConstant = -200, 
          centralGravity = 0.005,       
          springLength = 250,           
          springConstant = 0.05         
        ),
        stabilization = list(iterations = 150)
      ) %>%
      visEvents(
        stabilizationIterationsDone = "function() { this.setOptions({physics: false}); }",
        
        # Event 1: Click on a Node
        selectNode = "function(properties) { 
          Shiny.setInputValue('current_node_click', properties.nodes[0], {priority: 'event'}); 
        }",
        
        # Event 2: Click on Empty Space (Reset)
        click = "function(properties) { 
          if(properties.nodes.length === 0) {
            Shiny.setInputValue('empty_space_click', Math.random(), {priority: 'event'});
          }
        }"
      ) %>%
      visInteraction(
        hideNodesOnDrag = FALSE, hideEdgesOnDrag = FALSE, multiselect = TRUE
      )
  })
  
  # --- 4. Sync Graph Clicks to Dropdown ---
  observeEvent(input$current_node_click, {
    clicked_node <- input$current_node_click
    current_selection <- input$multi_node_select
    
    if (!(clicked_node %in% current_selection)) {
      updateSelectizeInput(session, "multi_node_select", 
                           selected = c(current_selection, clicked_node))
    }
  })
  
  # --- 5. Handle Empty Space Click (Reset) ---
  observeEvent(input$empty_space_click, {
    # Clear the dropdown. This will trigger the observer below to reset the graph style.
    updateSelectizeInput(session, "multi_node_select", selected = character(0))
  })
  
  # --- 6. Multi-Selection Highlighting Logic ---
  observeEvent(input$multi_node_select, {
    
    selected_ids <- input$multi_node_select
    dat <- graph_data()
    current_nodes <- dat$nodes
    current_edges <- dat$edges
    
    proxy <- visNetworkProxy("network_plot")
    
    # --- RESET LOGIC ---
    if (is.null(selected_ids) || length(selected_ids) == 0) {
      nodes_reset <- current_nodes %>% mutate(opacity = 1)
      edges_reset <- current_edges 
      
      visUpdateNodes(proxy, nodes = nodes_reset)
      visUpdateEdges(proxy, edges = edges_reset)
      
    } else {
      
      # --- HIGHLIGHT LOGIC ---
      incident_edges <- current_edges %>% 
        filter(from %in% selected_ids | to %in% selected_ids)
      
      neighbor_ids <- unique(c(incident_edges$from, incident_edges$to))
      nodes_to_keep_ids <- unique(c(selected_ids, neighbor_ids))
      edges_to_keep_ids <- incident_edges$id
      
      nodes_update <- current_nodes %>%
        mutate(
          color = ifelse(id %in% nodes_to_keep_ids, color, "rgba(200,200,200,0.3)"),
          opacity = ifelse(id %in% nodes_to_keep_ids, 1, 0.3)
        )
      
      edges_update <- current_edges %>%
        mutate(
          color = ifelse(id %in% edges_to_keep_ids, color, "rgba(220,220,220,0.1)"),
          width = ifelse(id %in% edges_to_keep_ids, width, 1)
        )
      
      visUpdateNodes(proxy, nodes = nodes_update)
      visUpdateEdges(proxy, edges = edges_update)
    }
  }, ignoreNULL = FALSE)
}

shinyApp(ui = ui, server = server)
