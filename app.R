library(shiny)
library(datamods)

# library(rsconnect)
# rsconnect::setAccountInfo(name='briankleiboeker', token='F659DF4BB6CE768CC7EB6E53272E64C9', secret='ECLyDsv90Prm15BEA/wLGJInXjSsgxU2JUaAQ/6m')
# deployApp(appDir = '/Users/briankleiboeker/Desktop/Lodhi lab/PexRAP_data_and_grant/20230105 3t3l1 kd pexrap 13c glucose lipidomics/shinyapp/structure_from_comp')

# ideas:
# 
# 3. find some way to include the analysis parameters in output csv?

import_copypaste_ui2 <- function (id, title = TRUE, name_field = TRUE) {
  ns <- NS(id)
  if (isTRUE(title)) {
    title <- tags$h4(i18n("Copy & paste data"), class = "datamods-title")
  }
  tags$div(class = "datamods-import", html_dependency_datamods(),
           title, tagAppendAttributes(textAreaInput(inputId = ns("data_pasted"),
                                                    label = i18n("Paste data here:"), height = "300px",
                                                    width = "100%", resize = "none",placeholder = "Data can be all 7 columns of data exactly as they come out of Thermo XCalibur (include column names):\n\nm/z	Int 	Rel.Int Theo.Mass	Delta	RDBE	Composition\n453.2815	341299.1	0.52	453.2816	-0.14	1.5	C20 [13]C H43 O7 N P \n453.3009	211562.2	0.32	453.2987	2.19 	0.5	C22 H46 O7 P \n454.2574	302408.3	0.46	454.2575	-0.14	1.5	C20 H41 O8 N P\n\nOr any other format with one column containing elemental compositions (minimum 2 columns, include column names still):\n\nm/z  	Composition\n453.2815	C20 [13]C H43 O7 N P \n453.3009	C22 H46 O7 P \n454.2574	C20 H41 O8 N P \n\nIf column names are not included, it will be assumed that the rightmost column contains compositions (minimum 2 columns)."), class = "shiny-input-container-inline"),
           if (isTRUE(name_field)) {
             textInput(inputId = ns("name"), label = NULL, placeholder = i18n("Enter sample name here"),
                       width = "100%")
           }, tags$div(id = ns("import-placeholder"), alert(id = ns("import-result"),
                                                            status = "info", 
                                                            tags$b(i18n("Nothing pasted yet!")),
                                                            i18n("Please copy and paste some data in the dialog box above."),
                                                            dismissible = TRUE)), uiOutput(outputId = ns("container_valid_btn"),
                                                                                           style = "margin-top: 20px;"))
}
environment(import_copypaste_ui2) <- asNamespace('datamods')
assignInNamespace("import_copypaste_ui", import_copypaste_ui2, ns = "datamods")

# define functions
extract_num_elements <- function(element_letter,composition){
  if(element_letter=="N"|element_letter=="n"){
    element_letter<-"N(?!a)"
  }
  if(element_letter=="C"|element_letter=="c"){
    element_letter<-"C(?!l)"
  }
  #if there is a square bracket (i.e. if it's an isotope), substitute in the double brackets so it works with regex below
  if(grepl("\\[",element_letter)){
    element_letter<-gsub("\\[","\\\\[",element_letter)
    element_letter<-gsub("\\]","\\\\]",element_letter)
  }
  
  num <- ifelse(!grepl(element_letter,composition,perl = T,ignore.case = T),
                #paste0(element_letter,"(?=\\s)|",element_letter,"(?=[0-9]+)|",element_letter,"$")
                0,  #return zero
                ifelse(!grepl(paste0("(?<=",element_letter,")[0-9]+|(?<=",element_letter,"\\s)[0-9]+"),composition,perl = T,ignore.case = T), #else, 
                       1,
                       as.numeric(regmatches(composition,regexpr(paste0("(?<=",element_letter,")[0-9]+|(?<=",element_letter,"\\s)[0-9]+"),composition,perl = T,ignore.case = T)))))
  return(num)
}
return_num_elements <- function(element_letter,column_of_df){
  sapply(column_of_df,function(x) extract_num_elements(element_letter,x))
}
assign_species <- function(c,h,o,n,p,na,cl,ion.mode,adducts,rdb,domain,lois,max.dbl.bnds){
  if(is.na(rdb) | is.null(adducts)){
    return(NA)
  }
  if(ion.mode == "neutral"){
    return(structure_from_exact_comp(c,h,o,n,p,c("nonacylglycerols","acylglycerols"),domain,lois,as.numeric(max.dbl.bnds)))
  }
  if(ion.mode == "pos.ion"){
    if("m.plus.sodium" %in% adducts){
      if(na > 0){
        return(structure_from_exact_comp(c,h,o,n,p,c("nonacylglycerols","acylglycerols"),domain,lois,as.numeric(max.dbl.bnds)))
        # look for sodium adduct species -- ALL
        # assign exact composition by changing nothing
        # regardless if find one or not, return the result here
      }
      # else, keep going
    }
    if("m.plus.ammonia" %in% adducts){
      if(n > 0){
        result <- structure_from_exact_comp(c,h-4,o,n-1,p,c("acylglycerols"),domain,lois,as.numeric(max.dbl.bnds))
        if(!is.na(result)){
          return(result)
        }
        # look for ammonia adduct species -- only DAG/TAG 
        # do this by subtracting 4 H's and 1 N and seeing if the species could exist
        # if find a match, return it, else, go look for m.plus.h ion if ("m.plus.h" %in% adducts)
        # but if we do find one here, add it to a list of things to return
      }
      # else, keep going
    }
    if("m.plus.h" %in% adducts){
      return(structure_from_exact_comp(c,h-1,o,n,p,c("nonacylglycerols"),domain,lois,as.numeric(max.dbl.bnds)))
      # look for lipid matches
      # do this by subtracting 1H and seeing if the species could exist
      # if it does exist -> return m.plus.h
      # if it doesn't exist -> return some clever paste0 showing what we searched for, beginning with 'Species not detected (traceback:'
    }
    return(NA)
    # if we make it here, return some clever paste0 showing what we searched for, beginning with 'Species not detected (traceback:'
    
  }
  if(ion.mode == "neg.ion"){
    if("m.plus.chloride" %in% adducts){
      if(cl > 0){
        return(structure_from_exact_comp(c,h,o,n,p,c("nonacylglycerols"),domain,lois,as.numeric(max.dbl.bnds)))
        # look for chloride adduct species
        # assign exact composition by changing nothing
        # regardless if find one or not, return the result here
      }
      # else, keep going
    }
    if("m.minus.h" %in% adducts){
      result <- structure_from_exact_comp(c,h+1,o,n,p,c("nonacylglycerols"),domain,lois,as.numeric(max.dbl.bnds))
      if(!is.na(result)){
        return(result)
      }
      # look for lipid matches
      # do this by adding 1H and seeing if the species could exist
      # if it does exist return m.minus.h
      # if it doesn't exist -> return some clever paste0 showing what we searched for, beginning with 'Species not detected (traceback:'
      
    }
    if("m.minus.2h" %in% adducts){
      
      result <- structure_from_exact_comp(c,h+2,o,n,p,c("nonacylglycerols"),domain,lois,as.numeric(max.dbl.bnds))
      if(!is.na(result) && grepl("CL",result)){
        return(result)
      }
      # look for lipid matches
      # do this by adding 1H and seeing if the species could exist
      # if it does exist return m.minus.h
      # if it doesn't exist -> return some clever paste0 showing what we searched for, beginning with 'Species not detected (traceback:'
      
    }
    return(NA)
  }
}
assign.el.vs.plasmalogen.vs.lyso <- function(c,o,rdb,max.o,general.class,base.carbons,lyso.carbon.cutoff,base.rdb){
  if(o == max.o){
    #not an EL or lyso- species
    return(paste0(general.class,"("))
  }else{
    #either an EL or lyso-species (but not both, as that would be o-2)
    if((c - base.carbons) <= lyso.carbon.cutoff){
      #lyso species
      return(paste0("L",general.class,"("))
    }else{
      #an EL
      if(rdb == (base.rdb- 1)){
        #a non-plasmalogen EL
        return(paste0(general.class,"(O-"))
      }else{
        #a plasmalogen
        return(paste0(general.class,"(P-"))
      }
      
    }
  }
}
structure_from_exact_comp <- function(c,h,o,n,p,which_to_look_for,domain,lois,max.dbl.bnds){
  if(domain == "bact"){
    lyso.carbon.cutoff.value <- 19
  }else{
    lyso.carbon.cutoff.value <- 27
  }
  rdb.equiv = ( c - (h/2) + ((n+p)/2) +1 )
  if(rdb.equiv %% 1 != 0 | (((2*(c+n))-h)/2) > 14 ){
    return(NA)
  }
  
  if( "acylglycerols" %in% which_to_look_for){
    if(n == 0 && p == 0 && o == 5 && c >= 35){
      #DAG
      if(rdb.equiv < 2 | !("dag" %in% lois)){
        return(NA)
      }

      n.double.bonds <-  (rdb.equiv - 2)
      n.fa.carbons <- (c - 3)
      
      if(n.double.bonds > max.dbl.bnds){
        return(NA)
      }
      
      return(paste0("DAG(",n.fa.carbons,":",n.double.bonds,")"))
      #DAG has base 2 RDB , base 3 carbons
    }
    if(n == 0 && p == 0 && o == 6 && c >= 49){
      #TAG
      if(rdb.equiv < 3 | !("tag" %in% lois)){
        return(NA)
      }
      
      n.double.bonds <-  (rdb.equiv - 3)
      n.fa.carbons <- (c - 3)
      
      if(n.double.bonds > max.dbl.bnds){
        return(NA)
      }
      
      return(paste0("TAG(",n.fa.carbons,":",n.double.bonds,")"))
      #TG has base 3 RDB , base 3 carbons
    }
    if(length(which_to_look_for)==1){
      return(NA)
    }
  }
  if("nonacylglycerols" %in% which_to_look_for){ #look for everything but TAG and DAG
    if(n == 0){
      if(p == 1){
        #pa, pg, or pi
        if(o >=7 && o<= 8){
          #PA has base 2 rdb , base 3 carbons
          if((o == 8 && rdb.equiv < 2) | (o == 7 && rdb.equiv < 1) | !("pa" %in% lois)){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 2)
          n.fa.carbons <- (c - 3)
          species <- assign.el.vs.plasmalogen.vs.lyso(c,o,rdb.equiv,8,"PA",3,lyso.carbon.cutoff.value,2)
          if(grepl("O-",species) | grepl("^L",species)){
            n.double.bonds <- n.double.bonds+1
          }
          
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          
          return(paste0(species,n.fa.carbons,":",n.double.bonds,")"))
          
        }
        if(o >=9 && o<= 10 ){
          #PG has base 2 rdb , base 6 carbons
          if((o == 10 && rdb.equiv < 2) | (o == 9 && rdb.equiv < 1) | !("pg" %in% lois)){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 2)
          n.fa.carbons <- (c - 6)

          species <- assign.el.vs.plasmalogen.vs.lyso(c,o,rdb.equiv,10,"PG",6,lyso.carbon.cutoff.value,2)
          if(grepl("O-",species) | grepl("^L",species)){
            n.double.bonds <- n.double.bonds+1
          }
          
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          
          return(paste0(species,n.fa.carbons,":",n.double.bonds,")"))
          
        }
        if(o >=12 && o<= 13 ){
          #PI has base 3 rdb , base 9 carbons
          if((o == 13 && rdb.equiv < 3) | (o == 12 && rdb.equiv < 2) | !("pi" %in% lois)){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 3)
          n.fa.carbons <- (c - 9)
          
          species <- assign.el.vs.plasmalogen.vs.lyso(c,o,rdb.equiv,13,"PI",9,lyso.carbon.cutoff.value,3)
          if(grepl("O-",species) | grepl("^L",species)){
            n.double.bonds <- n.double.bonds+1
          }
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          
          return(paste0(species,n.fa.carbons,":",n.double.bonds,")"))
          
        }
        return(NA)
      }
      if(c > 35 && p == 2){
        #some type of CL
        if(!("cl" %in% lois)){
          return(NA)
        }else{
        
        if(o == 15){
          #dilysoCL has 2 rdb, base 9 carbons
          if(rdb.equiv < 2){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 2)
          n.fa.carbons <- (c - 9)
          
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          
          return(paste0("dilysoCL(",n.fa.carbons,":",n.double.bonds,")"))
          
        }
        
        if(o == 16 && c > 50){
          #lysoCL has 3 rdb , base 9 carbons
          if(rdb.equiv < 3){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 3)
          n.fa.carbons <- (c - 9)
          
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          
          return(paste0("monolysoCL(",n.fa.carbons,":",n.double.bonds,")"))
          
        }
        
        if(o == 17 && c > 60){
          #CL has base 4 rdb , base 9 carbons
          if(rdb.equiv < 4){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 4)
          n.fa.carbons <- (c - 9)
          
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          
          return(paste0("CL(",n.fa.carbons,":",n.double.bonds,")"))
        }
        return(NA)
        }
      }
      if(p == 0 && o >= 2 && o <= 4 && c > 13){
        #some FFA probably
        if(o == 2 && c > 13 && c < 33){
          #plain old FFA
          if(rdb.equiv < 1 | !("ffa" %in% lois)){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 1)
          n.fa.carbons <- (c - 0)
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          return(paste0(n.fa.carbons,":",n.double.bonds,"-FA"))
        }
        
        if(o == 3 && c > 13 && c < 33){
          #hFFA
          if(rdb.equiv < 1 | !("ffa" %in% lois)){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 1)
          n.fa.carbons <- (c - 0)
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          return(paste0("h",n.fa.carbons,":",n.double.bonds,"-FA"))
        }
        
        if(o == 4 && c > 25 && c < 65){
          #FAHFA
          if(rdb.equiv < 2 | !("ffa" %in% lois)){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 2)
          n.fa.carbons <- (c - 0)
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          return(paste0(n.fa.carbons,":",n.double.bonds,"-FAHFA"))
        }
        return(NA)
      }
      return(NA)
    }
    if(n == 2){
      if(p == 1 && o >= 5 && o<= 7){
        #SM has 1 rdb, base 5 carbons
        if(rdb.equiv < 1 | !("sm" %in% lois)){
          return(NA)
        }
        n.double.bonds <-  (rdb.equiv - 1)
        n.fa.carbons <- (c - 5)
        if(n.double.bonds > max.dbl.bnds){
          return(NA)
        }
        return(paste0("SM(",n.fa.carbons,":",n.double.bonds,")"))
        
      }
      return(NA)
    }
    if(n == 1){
      if(p == 0 && o >= 3 && o <= 5){
        #if o == 4 it could be MAG+NH4 adduct
        #ceramides have 1 rdb, base 0 carbons
        if(rdb.equiv < 1 | !("cer" %in% lois)){
          return(NA)
        }
        n.double.bonds <-  (rdb.equiv - 1)
        n.fa.carbons <- (c)
        if(n.double.bonds > max.dbl.bnds){
          return(NA)
        }

        if(o == 3){
          if(n.double.bonds == 0){
            return(paste0("dihydroCer(",n.fa.carbons,":",n.double.bonds,")"))
          }else{
            return(paste0("Ceramide(",n.fa.carbons,":",n.double.bonds,")"))
          }
        }else{
          return(paste0("phytoCer(",n.fa.carbons,":",n.double.bonds,")"))
        }
      }
      if(p == 1){
        if(o >= 9 && o <= 10 ){
          #"PS has 3 rdb, base 6 carbons
          if((o == 10 && rdb.equiv < 3) | (o == 9 && rdb.equiv < 2) | !("ps" %in% lois)){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 3)
          n.fa.carbons <- (c - 6)
          species <- assign.el.vs.plasmalogen.vs.lyso(c,o,rdb.equiv,10,"PS",6,lyso.carbon.cutoff.value,3)
          if(grepl("O-",species) | grepl("^L",species)){
            n.double.bonds <- n.double.bonds+1
          }
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          return(paste0(species,n.fa.carbons,":",n.double.bonds,")"))
          
        }
        if(o>=7 && o<=8 ){
          if((c %% 2) == 0) {
            #PC has 2 rdb , base 8 carbons
            if((o == 8 && rdb.equiv < 2) | (o == 7 && rdb.equiv < 1)| !("pc" %in% lois)){
              return(NA)
            }
            n.double.bonds <-  (rdb.equiv - 2)
            n.fa.carbons <- (c - 8)
            species <- assign.el.vs.plasmalogen.vs.lyso(c,o,rdb.equiv,8,"PC",8,lyso.carbon.cutoff.value,2)
            if(grepl("O-",species) | grepl("^L",species)){
              n.double.bonds <- n.double.bonds+1
            }
            if(n.double.bonds > max.dbl.bnds){
              return(NA)
            }
            return(paste0(species,n.fa.carbons,":",n.double.bonds,")"))
            
          }else{
            #PE has 2 rdb , base 5 carbons
            if((o == 8 && rdb.equiv < 2) | (o == 7 && rdb.equiv < 1) | !("pe" %in% lois)){
              return(NA)
            }
            n.double.bonds <-  (rdb.equiv - 2)
            n.fa.carbons <- (c - 5)
            species <- assign.el.vs.plasmalogen.vs.lyso(c,o,rdb.equiv,8,"PE",5,lyso.carbon.cutoff.value,2)
            if(grepl("O-",species) | grepl("^L",species)){
              n.double.bonds <- n.double.bonds+1
            }
            if(n.double.bonds > max.dbl.bnds){
              return(NA)
            }
            return(paste0(species,n.fa.carbons,":",n.double.bonds,")"))
            
          }
        }
      }
      if(length(which_to_look_for)==1){
        return(NA)
      }
    }
  }
  return(NA)
}
most<-function(x){
  return( (length(which(x))/length(x)) > 0.95 )
}

# make lists
lipidclasslist <- c("DAG" = "dag","TAG" = "tag","PA" = "pa","PG" = "pg","PI" = "pi","CL" = "cl","SM" = "sm","Ceramide" = "cer","PC" = "pc","PE" = "pe","PS" = "ps","FFA" = "ffa")
negionclasses <- c("pe","cer","cl","pi","pg","pa","ps","ffa")
posionclasses <- c("pc","sm","tag","dag")
relevant.col.names <- c("mz","intensity","rel.intensity","theo.mass","delta","rdb.equiv","composition")

ui <- fluidPage(
  titlePanel("Assign lipid structure from composition"),
  tabsetPanel(
    id = "tabsid",
    tabPanel(value = 0,
             titlePanel("Enter a single composition"),
             sidebarLayout(
               sidebarPanel(  #inputs go here
                 #domain (eukaryote vs. prokaryote, as of now only affects lyso carbon cutoff variable)
                 selectInput("domain","Sample source",
                             c("Eukaryota" = "euk",
                               "Bacteria" = "bact"),
                             selected = "euk"),
                 #ionmode
                 selectInput("ionmode", "What mode were samples run in?",
                             c("negative ion mode" = "neg.ion",
                               "positive ion mode" = "pos.ion",
                               "N/A (search for structure as a neutral species)" = "neutral")),
                 conditionalPanel(
                   condition = "input.tabsid == 0  && input.ionmode == 'pos.ion'",
                   checkboxGroupInput("posionadduct", "Which ions/adducts would you like to look for?",
                                      c("[M+H]+" = "m.plus.h",
                                        "[M+NH4]+" = "m.plus.ammonia",
                                        "[M+Na]+" = "m.plus.sodium")),
                   selectInput("posionloi","Which lipid classes to search for (we advise against changing this)",
                               lipidclasslist,
                               selected = posionclasses,
                               multiple = T)
                 ),
                 conditionalPanel(
                   condition = "input.tabsid == 0  && input.ionmode == 'neg.ion'",
                   checkboxGroupInput("negionadduct", "Which ions/adducts would you like to look for?",
                                      c("[M-H]-" = "m.minus.h",
                                        "[M+Cl]-" = "m.plus.chloride",
                                        "[M-2H]2-" = "m.minus.2h")),
                   selectInput("negionloi","Which lipid classes to search for (we advise against changing this)",
                               lipidclasslist,
                               selected = negionclasses,
                               multiple = T)
                 ),
                 conditionalPanel(
                   condition = "input.tabsid == 0  && input.ionmode == 'neutral'",
                   selectInput("neutralloi","Which lipid classes to search for",
                               lipidclasslist,
                               selected = c(negionclasses,posionclasses),
                               multiple = T)
                 ),
                 sliderInput("maxdblbnds",
                             "Maximum number of acyl/alkyl double bonds allowed",
                             min = 0,max = 14,value = 14,step = 1,round = T,animate = F, ticks = F),
                 textInput("comp", "Enter composition", value = "", width = NULL, placeholder = "e.g. C36 H73 O8 N P"),
               ), 
               mainPanel( #outputs go here
                 htmlOutput("analyzed_data")
               )
             )
    ),
    
    tabPanel(value = 1,
             titlePanel("Paste all data for a sample"),
             fluidRow(
               column(
                 width = 6,
                 selectInput("bulkdomain","Sample source",
                             c("Eukaryota" = "euk",
                               "Bacteria" = "bact"),
                             selected = "euk"),
                 selectInput("bulkisotope", "Isotope being used",
                             c("None" = "none",
                               "[13]C" = "carbon_isotope",
                               "[2]H" = "hydrogen_isotope")),
                 selectInput("bulkionmode", "What mode were samples run in?",
                             c("negative ion mode" = "neg.ion",
                               "positive ion mode" = "pos.ion",
                               "N/A (search for structure as a neutral species)" = "neutral")),
                 conditionalPanel(
                   condition = "input.tabsid == 1  && input.bulkionmode == 'pos.ion'",
                   checkboxGroupInput("bulkposionadduct", "Which ions/adducts would you like to look for?",
                                      c("[M+H]+" = "m.plus.h",
                                        "[M+NH4]+" = "m.plus.ammonia",
                                        "[M+Na]+" = "m.plus.sodium")),
                   selectInput("bulkposionloi","Which lipid classes to search for (we advise against changing this)",
                               lipidclasslist,
                               selected = posionclasses,
                               multiple = T)
                 ),
                 conditionalPanel(
                   condition = "input.tabsid == 1  && input.bulkionmode == 'neg.ion'",
                   checkboxGroupInput("bulknegionadduct", "Which ions/adducts would you like to look for?",
                                      c("[M-H]-" = "m.minus.h",
                                        "[M+Cl]-" = "m.plus.chloride",
                                        "[M-2H]2-" = "m.minus.2h")),
                   selectInput("bulknegionloi","Which lipid classes to search for (we advise against changing this)",
                               lipidclasslist,
                               selected = negionclasses,
                               multiple = T)
                 ),
                 conditionalPanel(
                   condition = "input.tabsid == 1  && input.bulkionmode == 'neutral'",
                   selectInput("bulkneutralloi","Which lipid classes to search for",
                               lipidclasslist,
                               selected = c(negionclasses,posionclasses),
                               multiple = T)
                 ),
                 radioButtons("returnAll", "Would you like to...", 
                              choiceNames = c("Return all peaks",
                                              "Only return peaks with a species assigned"),
                              choiceValues = c("returnallpeaks","returnassignedpeaks"),
                              inline = TRUE, selected = ("returnallpeaks")),
                 numericInput("deltacutoff", "abs(Delta) cutoff:", 10, min = 0.1, max = 1000,value = 1000),
                 sliderInput("bulkmaxdblbnds",
                             "Maximum number of acyl/alkyl double bonds allowed",
                             min = 0,max = 14,value = 14,step = 1,round = T,animate = F, ticks = F),
                 import_copypaste_ui2("myid",
                                      title = ""
                                      )
               ),
               column(
                 width = 6,
                 tags$b("Import status:"),
                 htmlOutput(outputId = "status"),
                 #verbatimTextOutput(outputId = "status"),
                 tags$b("Sample name:"),
                 verbatimTextOutput(outputId = "name"),
                 tags$b("Verify that data was input & read correctly (just make sure columns align):"),
                 verbatimTextOutput(outputId = "data"),
                 downloadButton("download", "Download results")
               )
             )
    )
  )
)

server <- function(input, output, session) {
  ############
  observeEvent(input$ionmode,
               {if(input$ionmode == "pos.ion"){
                 #reset negative and neutral
                 updateCheckboxGroupInput(session = session,
                                          "negionadduct",
                                          "Which ions/adducts would you like to look for?",
                                          c("[M-H]-" = "m.minus.h",
                                            "[M+Cl]-" = "m.plus.chloride",
                                            "[M-2H]2-" = "m.minus.2h"),
                                          selected=character(0))
                 updateSelectInput(session = session, "negionloi",
                                   selected = negionclasses)
                 updateSelectInput(session = session, "neutralloi",
                                   selected = c(negionclasses,posionclasses))
               }
               if(input$ionmode == "neg.ion"){
                 updateCheckboxGroupInput(session = session,
                                          "posionadduct",
                                          "Which ions/adducts would you like to look for?",
                                          c("[M+H]+" = "m.plus.h",
                                            "[M+NH4]+" = "m.plus.ammonia",
                                            "[M+Na]+" = "m.plus.sodium"),
                                          selected=character(0))
                 updateSelectInput(session = session, "posionloi",
                                   selected = posionclasses)
                 updateSelectInput(session = session, "neutralloi",
                                   selected = c(negionclasses,posionclasses))
               }
               if(input$ionmode == "neutral"){
                  updateCheckboxGroupInput(session = session,
                                           "negionadduct",
                                           "Which ions/adducts would you like to look for?",
                                           c("[M-H]-" = "m.minus.h",
                                             "[M+Cl]-" = "m.plus.chloride",
                                             "[M-2H]2-" = "m.minus.2h"),
                                           selected=character(0))
                  updateSelectInput(session = session, "negionloi",
                                    selected = negionclasses)
                  updateCheckboxGroupInput(session = session,
                                           "posionadduct",
                                           "Which ions/adducts would you like to look for?",
                                           c("[M+H]+" = "m.plus.h",
                                             "[M+NH4]+" = "m.plus.ammonia",
                                             "[M+Na]+" = "m.plus.sodium"),
                                           selected=character(0))
                  updateSelectInput(session = session, "posionloi",
                                    selected = posionclasses)
                }
                 
                },
               ignoreInit = T)
  observeEvent(input$bulkionmode,
               {if(input$bulkionmode == "pos.ion"){
                 updateCheckboxGroupInput(session = session,
                                          "bulknegionadduct",
                                          "Which ions/adducts would you like to look for?",
                                          c("[M-H]-" = "m.minus.h",
                                            "[M+Cl]-" = "m.plus.chloride",
                                            "[M-2H]2-" = "m.minus.2h"),
                                          selected=character(0))
                 updateSelectInput(session = session, "bulknegionloi",
                                   selected = negionclasses)
                 updateSelectInput(session = session, "bulkneutralloi",
                                   selected = c(negionclasses,posionclasses))
                 
               }
                 if(input$bulkionmode == "neg.ion"){
                 updateCheckboxGroupInput(session = session,
                                          "bulkposionadduct",
                                          "Which ions/adducts would you like to look for?",
                                          c("[M+H]+" = "m.plus.h",
                                            "[M+NH4]+" = "m.plus.ammonia",
                                            "[M+Na]+" = "m.plus.sodium"),
                                          selected=character(0))
                 updateSelectInput(session = session, "bulkposionloi",
                                   selected = posionclasses)
                 updateSelectInput(session = session, "bulkneutralloi",
                                   selected = c(negionclasses,posionclasses))
               }
                 if(input$bulkionmode == "neutral"){
                   updateCheckboxGroupInput(session = session,
                                            "bulknegionadduct",
                                            "Which ions/adducts would you like to look for?",
                                            c("[M-H]-" = "m.minus.h",
                                              "[M+Cl]-" = "m.plus.chloride",
                                              "[M-2H]2-" = "m.minus.2h"),
                                            selected=character(0))
                   updateSelectInput(session = session, "bulknegionloi",
                                     selected = negionclasses)
                   updateCheckboxGroupInput(session = session,
                                            "bulkposionadduct",
                                            "Which ions/adducts would you like to look for?",
                                            c("[M+H]+" = "m.plus.h",
                                              "[M+NH4]+" = "m.plus.ammonia",
                                              "[M+Na]+" = "m.plus.sodium"),
                                            selected=character(0))
                   updateSelectInput(session = session, "bulkposionloi",
                                     selected = posionclasses)
                 }
                 },
               ignoreInit = T)
  
  output$analyzed_data <- renderUI({
    c <- extract_num_elements("C",input$comp)
    h <- extract_num_elements("H",input$comp)
    o <- extract_num_elements("O",input$comp)
    n <- extract_num_elements("N",input$comp)
    p <- extract_num_elements("P",input$comp)
    na <- extract_num_elements("Na",input$comp)
    cl <- extract_num_elements("Cl",input$comp)
    
    str1 <- paste0("Carbon: ",c)
    str2 <- paste0("Hydrogen: ",h)
    str3 <- paste0("Oxygen: ",o)
    str4 <- paste0("Nitrogen: ",n)
    str5 <- paste0("Phosphorus: ",p)
    str6 <- paste0("Sodium: ",na)
    str7 <- paste0("Chlorine: ",cl)
    
    str8<-""
    str9<-""
    rdb.equiv = ( c - (h/2) + ((n+p)/2) +1 )
    str10<-paste0("RDB equiv. = ", rdb.equiv)
    str11<-""

    if(input$ionmode == "pos.ion"){
      adducts <- input$posionadduct
      lois <- input$posionloi
    }
    if(input$ionmode == "neg.ion"){
      adducts <- input$negionadduct
      lois <- input$negionloi
    }
    if(input$ionmode == "neutral"){
      adducts <- c("m.plus.h","m.plus.ammonia","m.plus.sodium","m.minus.h","m.plus.chloride","m.minus.2h")
      lois <- input$neutralloi
    }
    
    species.result <- assign_species(c,h,o,n,p,na,cl,input$ionmode,adducts,rdb.equiv,input$domain,lois,input$maxdblbnds)
    
    if(!is.na(species.result)){
      str12<-paste0("Species is: ","<b>",species.result)
    }else{
      str12<-paste0("Species not detected")
    }
    
    if(input$comp == ""){
      url <- a("this colab notebook.", href="https://colab.research.google.com/drive/1rQz5ivnjUGu7ui5jsQAHEHZbAynC3xRP?usp=sharing", target="_blank")
      tagList(
        p(HTML(paste("","This webapp was designed for use with high-resolution Orbitrap ESI-MS data with the goal of accurately assigining lipid class and structure to user-provided elemental compositions, but should also be compatible with most any mass spectromerty data so long as the elemental compositions and column names are provided in a standard format.","","Report problems, suggestions, or other feedback to bkleiboeker [at] wustl.edu","", sep = '<br/>'))),
        p("Need to assign elemental compositions to m/z values after realigning m/z's? Check out ",url),
        p(HTML(paste("","","Input an elemental composition to begin",sep = '<br/>')))
      )
    }else{
      HTML(paste("Report problems, suggestions, or other feedback to bkleiboeker [at] wustl.edu"," ",str1, str2, str3, str4, str5, str6, str7,str8,str9,str10,str11,str12,sep = '<br/>'))
    }
    
  })
  
  ############
  imported <- import_copypaste_server("myid",
                                      btn_show_data = F,
                                      trigger_return = "change")
  output$status <- renderPrint({
    if(!is.null(imported$status())){
      if(imported$status() == "success"){
        if(!("composition" %in% colnames(data.full()))){
          if(length(which(colnames(data.full()) %in% relevant.col.names)) == 0 ){
            HTML(paste0("Import successful, but no column names detected. ", "Assuming that the rightmost non-empty column contains compositions. Filter output by abs(delta) cutoff feature is not available.",sep = '<br/>'))
          }else{
            HTML(paste0("Import successful, but composition column not detected. ", "Assuming that the rightmost non-empty column contains compositions. ",ifelse("delta" %in% colnames(data.full()),"Delta column detected, so filtering by abs(delta) cutoff feature is enabled. ","Delta column not detected, so filtering by abs(delta) cutoff feature is disabled. "),sep = '<br/>'))
          }
          
        }else{
          HTML(paste0("Import successful. Detected columns: ",paste( colnames(data.full())[colnames(data.full()) != "structure"],collapse = ', '),ifelse("delta" %in% colnames(data.full()),". Delta column detected, so filtering by abs(delta) cutoff feature is enabled. ",". Delta column not detected, so filtering by abs(delta) cutoff feature is disabled. ")))
        }
      }else{
        imported$status()
      }
    }else{
      imported$status()
    }

  })
  output$name <- renderPrint({
    imported$name()
  })
 
  ############
  data.full <- reactive({
    if(!is.null(imported$data())){
      df<-imported$data()
      
      while( most(is.na(df[,ncol(df)])) | most(df[,ncol(df)] == "") | most(is.null(df[,ncol(df)]))){
        df<-df[,1:(ncol(df)-1)]
      }
      while( most(is.na(df[,1])) | most(df[,1] == "") | most(is.null(df[,1]))){
        df<-df[,2:ncol(df)]
      }
      
      colnames(df)[(unlist(lapply(colnames(df),grepl,pattern = "m",ignore.case = T,perl = T)) & unlist(lapply(colnames(df),grepl,pattern = "z",ignore.case = T,perl = T)))]<-"mz"
      colnames(df)[(unlist(lapply(colnames(df),grepl,pattern = "int",ignore.case = T,perl = T)) & !unlist(lapply(colnames(df),grepl,pattern = "rel",ignore.case = T,perl = T)))]<-"intensity"
      colnames(df)[unlist(lapply(colnames(df),grepl,pattern = "rel",ignore.case = T,perl = T))]<-"rel.intensity"
      colnames(df)[unlist(lapply(colnames(df),grepl,pattern = "theo",ignore.case = T,perl = T))]<-"theo.mass"
      colnames(df)[( unlist(lapply(colnames(df),grepl,pattern = "delta",ignore.case = T,perl = T)) | unlist(lapply(colnames(df),grepl,pattern = "mmu",ignore.case = T,perl = T)) | unlist(lapply(colnames(df),grepl,pattern = "amu",ignore.case = T,perl = T)) | unlist(lapply(colnames(df),grepl,pattern = "ppm",ignore.case = T,perl = T)) )]<-"delta"
      colnames(df)[unlist(lapply(colnames(df),grepl,pattern = "rdb",ignore.case = T,perl = T))]<-"rdb.equiv"
      colnames(df)[unlist(lapply(colnames(df),grepl,pattern = "comp",ignore.case = T,perl = T))]<-"composition"
      
      if("composition" %in% colnames(df)){
         df$c <- return_num_elements("C",df$composition)
         df$h <- return_num_elements("H",df$composition)
         df$o <- return_num_elements("O",df$composition)
         df$n <- return_num_elements("N",df$composition)
         df$p <- return_num_elements("P",df$composition)
         df$na <- return_num_elements("Na",df$composition)
         df$cl <- return_num_elements("Cl",df$composition)
         if(input$bulkisotope != "none"){
           if(input$bulkisotope == "carbon_isotope"){
             df$c <- df$c + return_num_elements("[13]C",df$composition)
           }else{
             df$h <- df$h + return_num_elements("[2]H",df$composition)
           }
         }
      }else{
        ncols <- ncol(df)
        df$c <- return_num_elements("C",df[,ncols])
        df$h <- return_num_elements("H",df[,ncols])
        df$o <- return_num_elements("O",df[,ncols])
        df$n <- return_num_elements("N",df[,ncols])
        df$p <- return_num_elements("P",df[,ncols])
        df$na <- return_num_elements("Na",df[,ncols])
        df$cl <- return_num_elements("Cl",df[,ncols])
        if(input$bulkisotope != "none"){
          if(input$bulkisotope == "carbon_isotope"){
            df$c <- df$c + return_num_elements("[13]C",df[,ncols])
          }else{
            df$h <- df$h + return_num_elements("[2]H",df[,ncols])
          }
        }
      }
 
      if(input$bulkionmode == "pos.ion"){
        bulkadducts <- input$bulkposionadduct
        bulklois <- input$bulkposionloi
      }
      if(input$bulkionmode == "neg.ion"){
        bulkadducts <- input$bulknegionadduct
        bulklois <- input$bulknegionloi
      }
      if(input$bulkionmode == "neutral"){
        bulkadducts <- c("m.plus.h","m.plus.ammonia","m.plus.sodium","m.minus.h","m.plus.chloride","m.minus.2h")
        bulklois <- input$bulkneutralloi
      }
      
      df$structure <- mapply(assign_species,df$c,df$h,df$o,df$n,df$p,df$na,df$cl,input$bulkionmode,rep(list(bulkadducts),nrow(df)),( df$c - (df$h/2) + ((df$n+df$p)/2) +1 ),input$bulkdomain,rep(list(bulklois),nrow(df)),input$bulkmaxdblbnds,SIMPLIFY = T)
      df$structure <- ifelse(is.na(df$structure),"",df$structure)
      df<-df[, !(colnames(df) %in% c("c","h","o","n","p","na","cl"))]
      df
      
    }else{
      imported$data()
    }
    
  })
  

  output$data <- renderPrint({
    if(!is.null(imported$data()) ){
      if("delta" %in% colnames(data.full())){
        if(input$returnAll == "returnassignedpeaks"){
          data.full()[ifelse( !is.na(data.full()$delta) & abs(data.full()$delta) <= input$deltacutoff & data.full()$structure != "" ,T,F) ,] |> head(n=10)
        }else{
          data.full()[ifelse( !is.na(data.full()$delta) & abs(data.full()$delta) <= input$deltacutoff ,T,F) ,] |> head(n=10)
        }
      }else{
        if(input$returnAll == "returnassignedpeaks"){
          data.full()[ifelse(data.full()$structure != "" ,T,F) ,] |> head(n=10)
        }else{
          data.full() |> head(n=10)
        }
      }
    }else{
      data.full()
    }
  })
  
  output$download <- downloadHandler(
    filename = function() {
      paste0(imported$name(),".csv")
    },
    content = function(file) {
      if(!is.null(imported$data()) ){
        if("delta" %in% colnames(data.full())){
          if(input$returnAll == "returnassignedpeaks"){
            write.csv(data.full()[ifelse( !is.na(data.full()$delta) & abs(data.full()$delta) <= input$deltacutoff & data.full()$structure != "" ,T,F) ,],file,row.names = F)
          }else{
            write.csv(data.full()[ifelse( !is.na(data.full()$delta) & abs(data.full()$delta) <= input$deltacutoff ,T,F) ,],file,row.names = F)
          }
        }else{
          if(input$returnAll == "returnassignedpeaks"){
            write.csv(data.full()[ifelse(data.full()$structure != "" ,T,F) ,],file,row.names = F) 
          }else{
            write.csv(data.full() ,file,row.names = F)
          }
        }
      }else{
        write.csv( data.full(),file,row.names = F)
      }

    }
  )
  
}

shinyApp(ui = ui, server = server)
