# Definição da UI
ui <- tagList(
  
  # ---- Google Analytics 4 (GTAG) ----
  tags$head(
    tags$title("Multi-omic Oncometabolism GPS Shiny"),
    HTML("
      <!-- Google tag (gtag.js) -->
      <script async src='https://www.googletagmanager.com/gtag/js?id=G-0T0KK6YE0T'></script>
      <script>
        window.dataLayer = window.dataLayer || [];
        function gtag(){ dataLayer.push(arguments); }
        gtag('js', new Date());
        // Desliga page_view automático (SPA Shiny)
        gtag('config', 'G-0T0KK6YE0T', { 'send_page_view': false, 'anonymize_ip': true });

        // Função para enviar page_view manual
        function sendPageView(titleSuffix) {
          var loc = window.location.href;
          var title = titleSuffix ? ('Oncometabolism | ' + titleSuffix) : 'Oncometabolism';
          gtag('event', 'page_view', {
            page_title: title,
            page_location: loc
          });
        }

        // Dispara um page_view na carga inicial
        document.addEventListener('DOMContentLoaded', function() {
          sendPageView();
        });

        // Dispara page_view a cada troca de aba no navbar do Shiny
        document.addEventListener('DOMContentLoaded', function() {
          $(document).on('shown.bs.tab', 'a[data-toggle=\"tab\"], a[data-bs-toggle=\"tab\"]', function (e) {
            var tabText = $(e.target).text().trim();
            sendPageView(tabText);
          });
        });
      </script>
      <noscript>Google Analytics requires JavaScript enabled.</noscript>
    ")
  ),
  
  # ----------------- Tela de carregamento -----------------
  div(
    id = "loading_screen",
    style = "position: fixed; 
             top: 0; left: 0; 
             width: 100%; height: 100%; 
             background-color: white; 
             z-index: 9999; 
             display: flex; 
             flex-direction: column; 
             justify-content: center; 
             align-items: center;",
    h3("Loading Multi-omic Oncometabolism GPS...", style = "color: #386257; font-size: 24px;"),
    tags$style("
      @keyframes progress {
        0% { width: 0%; }
        100% { width: 100%; }
      }
    "),
    div(
      style = "width: 50%; background-color: #ddd; height: 10px; border-radius: 5px; overflow: hidden;",
      div(style = "height: 100%; width: 0%; background-color: #386257; animation: progress 30s linear;")
    )
  ),
  
  # Script para remover a tela de carregamento após 30 segundos
  tags$script("
    setTimeout(function() {
      document.getElementById('loading_screen').style.display = 'none';
    }, 30000);
  "),
  
  # ----------------- Interface principal -----------------
  navbarPage(
    title = div("Multi-omic Oncometabolism GPS Shiny"),
    inverse = TRUE,
    
    # ---------- CSS ----------
    header = tags$style(HTML("
      .navbar { background-color: #386257; }
      .navbar .navbar-brand { color: white; font-size: 18px; }
      .navbar .navbar-nav > li > a { color: white; }
      .navbar .navbar-nav > li > a:hover { background-color: #233D36; }
      .navbar .navbar-nav .active > a { background-color: #1D332D; }
      /* Centraliza verticalmente o ícone da aba do GitHub */
      .navbar .navbar-nav > li > a[data-value='github_link'] {
        padding-top: 12px;   /* ajuste fino da altura */
        padding-bottom: 12px;
      }
      /* Tamanho e alinhamento do ícone */
      .github-icon {
        font-size: 22px;     /* ajuste o tamanho que preferir */
        line-height: 1; 
        vertical-align: middle;
      }
    ")),
    
    # -------- Abas ----------
    tabPanel(
      "Home",
      fluidRow(
        column(
          12, align = "center",
          h3("Welcome to Multi-omic Oncometabolism GPS Shiny", style = "font-size: 32px;"),
          p(
            "Explore results from Omic layer- and metabolic pathway-specific signatures for pan-cancer biomarker and therapeutic target discovery.",
            style = "font-size: 18px;"
          ),
          img(src = "Figure_1.png", style = "max-width: 100%; height: 100%;")
        )
      ),
      br(), br(), hr(),
      fluidRow(
        column(
          12, align = "center",
          h4("Visitor Count", style = "color: #386257; font-size: 20px;"),
          textOutput("visitor_count")
        )
      ),
      br()
    ),
    
    tabPanel("User Manual", mod_user_manual_ui("manual")),
    
    # tabPanel("AI Podcast and Presentation",
    #          fluidRow(
    #            column(12, align = "center",
    #                   br(), br(),
    #                   
    #                   # Podcast title with headphones icon
    #                   tags$h3(tags$i(class = "fa fa-headphones"), " Podcast", 
    #                           style = "font-size: 50px; font-weight: bold;"),
    #                   
    #                   # Playable audio element (Fixed File Path)
    #                   tags$audio(controls = TRUE, 
    #                              tags$source(src = "Podcast.m4a", type = "audio/wav")),
    #                   
    #                   br(), br(),
    #                   
    #                   # Styled Podcast Button (Fixed File Path)
    #                   tags$a(href = "Podcast.m4a", 
    #                          download = "Podcast.m4a", 
    #                          class = "btn btn-lg",
    #                          style = "background-color: #2c3e50; color: white; border-color: #2c3e50; padding: 10px 20px; font-size: 18px; border-radius: 5px;",
    #                          "Download our Podcast")
    #            )
    #          )
    # ),
    
    tabPanel("AI Podcast and Presentation",
             fluidRow(
               column(12, align = "center",
                      br(), br(),
                      
                      # Podcast title with headphones icon
                      tags$h3(tags$i(class = "fa fa-headphones"), " Podcast", 
                              style = "font-size: 50px; font-weight: bold;"),
                      
                      # Playable audio element
                      tags$audio(
                        controls = TRUE, 
                        tags$source(src = "Podcast.m4a", type = "audio/mp4")
                      ),
                      
                      br(), br(),
                      
                      # Styled Podcast Button
                      tags$a(href = "Podcast.m4a", 
                             download = "Podcast.m4a", 
                             class = "btn btn-lg",
                             style = "background-color: #2c3e50; color: white; border-color: #2c3e50; padding: 10px 20px; font-size: 18px; border-radius: 5px;",
                             "Download our Podcast"),
                      
                      br(), br(), br(),
                      
                      # VIDEO SECTION
                      tags$h3(tags$i(class = "fa fa-video-camera"), " Presentation Video",
                              style = "font-size: 45px; font-weight: bold;"),
                      
                      tags$video(
                        width = "800px",
                        controls = NA,
                        tags$source(src = "Presentation.mp4", type = "video/mp4")
                      )
               )
             )
    ),
    
    
    # tabPanel("Exploration of molecular targets", mod_atlas_ui("atlas")),
    
    tabPanel("Search Your Gene", mod_search_your_target_ui("search_your_target")),
    
    navbarMenu(
      "Filter Signatures",
      tabPanel("CNV-Specific Signatures",         mod_signature_ui("cnv_signature")),
      tabPanel("mRNA-Specific Signatures",        mod_signature_ui("gene_signature")),
      tabPanel("Methylation-Specific Signatures", mod_signature_ui("methylation_signature")),
      tabPanel("Mutation-Specific Signatures",    mod_signature_ui("mutation_signature")),
      tabPanel("miRNA-Specific Signatures",       mod_signature_ui("mirna_signature")),
      tabPanel("Protein-Specific Signatures",     mod_signature_ui("protein_signature")),
      tabPanel("Transcript-Specific Signatures",  mod_signature_ui("transcript_signature"))
    ),
    
    # tabPanel(
    #   "Top Signatures",
    #   mod_the_best_signature_ui("top_ranked_signatures")
    # ),
    
    # tabPanel(
    #   "Meaningful Interaction",
    #   mod_meaningful_interaction_ui("Mean_int")
    # ),
    # 
    # tabPanel(
    #   "Regulatory Circuitries",
    #   mod_regulatory_circuitry_ui("reg_circ")
    # ),
    
    tabPanel(
      title = "Regulatory Circuitries",
      mod_interaction_network_ui("interaction_net")
    ),
    
    # # Example inside a navbarPage/tabsetPanel
    # tabPanel(
    #   "Custom Signatures",
    #   mod_custom_signatures_ui("custom_signatures")
    # ),
    
    navbarMenu(
      "Analysis and Plotting",
      tabPanel("Correlation Analysis",        mod_correlation_analysis_ui("correlation_analysis")),
      tabPanel("Tumor vs Normal Analysis",    mod_tumor_normal_analysis_ui("tumor_normal_analysis")),
      tabPanel("Cox Analysis",                mod_cox_analysis_ui("cox_analysis")),
      tabPanel("Survival Analysis",           mod_survival_analysis_ui("survival_analysis")),
      tabPanel("Immune Infiltrates Analysis", mod_infiltrates_analysis_ui("infiltrates_analysis"))
    ),
    
    tabPanel("About us", mod_developers_ui("Developers")),
    
    tabPanel(
      title = "Citation",
      value = "citation_tab",
      fluidPage(
        h3("How to Cite"),
        
        p("If you use OncoMetabolismGPS or results generated by this application, please cite:"),
        
        tags$div(
          style = "margin-left: 10px;",
          p(
            "Nogueira, H. A. C.; Souza, E. R.; Lopes, V. S.; Medina-Acosta, E. (2025) ",
            em("A Multi-Omic Atlas of Convergent and Divergent Metabolic Regulatory Circuitries in Cancer."),
            " Preprint. DOI: ",
            a(
              "10.1101/2025.11.15.688631",
              href = "https://doi.org/10.1101/2025.11.15.688631",
              target = "_blank"
            )
          )
        ),
        
        h4("BibTeX"),
        tags$pre(
          "@article{Nogueira2025MultiOmicAtlas,
  title={A Multi-Omic Atlas of Convergent and Divergent Metabolic Regulatory Circuitries in Cancer},
  author={Nogueira, Higor A. C. and Souza, E. R. and Lopes, V. S. and Medina-Acosta, E.},
  year={2025},
  journal={Preprint},
  doi={10.1101/2025.11.15.688631}
}"
        )
      )
    ),
    
    
    # -------- ABA "GITHUB" (sem <a> interno) --------
    tabPanel(
      title = HTML('<span class="fa fa-github github-icon"></span>'),  # FA4
      value = "github_link"
    )
  ),
  
  # ---------- JS: clique na aba abre o GitHub em nova aba ----------
  tags$script(HTML("
    $(document).on('click', 'a[data-value=\"github_link\"]', function(e){
      e.preventDefault();
      window.open('https://github.com/BioCancerInformatics/Multi-omic-Oncometabolism-GPS', '_blank');
    });
  "))
)
