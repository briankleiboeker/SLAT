library(shiny)
library(datamods)
library(DT)
library(ggplot2)
library(rclipboard)
library(plotly)

#todo: think of ways for user to add new lipid species (temp or permanent)

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
import_copypaste_ui3 <- function (id, title = TRUE, name_field = TRUE) {
  ns <- NS(id)
  if (isTRUE(title)) {
    title <- tags$h4(i18n("Copy & paste data"), class = "datamods-title")
  }
  tags$div(class = "datamods-import", html_dependency_datamods(),
           title, tagAppendAttributes(textAreaInput(inputId = ns("data_pasted"),
                                                    label = i18n("Paste data here:"), height = "300px",
                                                    width = "100%", resize = "none",placeholder = "Data must contain at least 2 columns (more is ok). Please include column names:\n\nm/z	         Intensity \n453.2815	 341299.1 \n453.3009	 211562.2 \n454.2574 302408.3\n\nThe column containing m/z values should be named 'mz' or 'm/z' (case insensitive). If column names are not included, it will be assumed that the leftmost column contains m/z values (minimum 2 columns still)."), class = "shiny-input-container-inline"),
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
environment(import_copypaste_ui3) <- asNamespace('datamods')
assignInNamespace("import_copypaste_ui", import_copypaste_ui3, ns = "datamods")

# define functions
extract_num_elements <- function(element_letter,composition){
  if(element_letter=="N"|element_letter=="n"){
    element_letter<-"(?<!I|M|Z|S|R)N(?!a|s|d|e|p|i|b|o)"
  }
  if(element_letter=="C"|element_letter=="c"){
    element_letter<-"(?<!A|T|S)C(?!a|d|f|e|s|l|r|o|u|m)"
  }
  if(element_letter=="H"|element_letter=="h"){
    element_letter<-"(?<!T|R)H(?!g|f|s|e|o)"
  }
  if(element_letter=="O"|element_letter=="o"){
    element_letter<-"(?<!H|C|P|N)O(?!s)"
  }
  if(element_letter=="P"|element_letter=="p"){
    element_letter<-"(?<!N)P(?!b|d|t|u|o|r|m|a)"
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
        return(structure_from_exact_comp(c,h,o,n,p,c("nonacylglycerols","acylglycerols","galactosyldiacylglycerols"),domain,lois,as.numeric(max.dbl.bnds)))
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
      # regardless if find one or not, return the result here
    }
    return(NA)
    # if we make it here, return some clever paste0 showing what we searched for, beginning with 'Species not detected (traceback:'
    
  }
  if(ion.mode == "neg.ion"){
    if("m.plus.chloride" %in% adducts){
      if(cl > 0){
        return(structure_from_exact_comp(c,h,o,n,p,c("nonacylglycerols","galactosyldiacylglycerols"),domain,lois,as.numeric(max.dbl.bnds)))
        # look for chloride adduct species
        # assign exact composition by changing nothing
        # regardless if find one or not, return the result here
      }
      # else, keep going
    }
    if("m.minus.h" %in% adducts){
      result <- structure_from_exact_comp(c,h+1,o,n,p,c("nonacylglycerols","galactosyldiacylglycerols"),domain,lois,as.numeric(max.dbl.bnds))
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
        #this is true no matter what
        return(paste0(general.class,"(O-"))
      }else{
        # there are one or more double bonds -- these can be adjacent to ether linkage (plasmalogen) or anywhere else
        if(general.class %in% c("PC","PI")){
          return(paste0(general.class,"(O-"))
        }else{
          return(paste0(general.class,"(P-"))
        }
        
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
  
  if(rdb.equiv %% 1 != 0 ){
    return(NA)
  }
  
  #put an if(strain == "bacteroidetes") here
  #look for 
  
  if( "acylglycerols" %in% which_to_look_for){
    if(n == 0 && p == 0 && o == 5 && c >= 35){
      #DAG
      #DAG has base 2 RDB , base 3 carbons
      if(rdb.equiv < 2 | !("dag" %in% lois)){
        return(NA)
      }
      n.double.bonds <-  (rdb.equiv - 2)
      n.fa.carbons <- (c - 3)
      if(n.double.bonds > max.dbl.bnds){
        return(NA)
      }
      return(paste0("DAG(",n.fa.carbons,":",n.double.bonds,")"))
    }
    if(n == 0 && p == 0 && o == 6 && c >= 49){
      #TAG
      #TAG has base 3 RDB , base 3 carbons
      if(rdb.equiv < 3 | !("tag" %in% lois)){
        return(NA)
      }
      n.double.bonds <-  (rdb.equiv - 3)
      n.fa.carbons <- (c - 3)
      if(n.double.bonds > max.dbl.bnds){
        return(NA)
      }
      return(paste0("TAG(",n.fa.carbons,":",n.double.bonds,")"))
    }
    if(length(which_to_look_for)==1){
      return(NA)
    }
  }
  if("galactosyldiacylglycerols" %in% which_to_look_for){
    if(n == 0 && p == 0 && c > 29 && (o == 10 | o == 15)){
      if(o == 15){
        if(rdb.equiv < 4 | !("gdg" %in% lois)){
          return(NA)
        }
        n.double.bonds <-  (rdb.equiv - 4)
        n.fa.carbons <- (c - 15)
        if(n.double.bonds > max.dbl.bnds){
          return(NA)
        }
        return(paste0("DGDG(",n.fa.carbons,":",n.double.bonds,")"))
      }
    if(o == 10){
      if(rdb.equiv < 3 | !("gdg" %in% lois)){
        return(NA)
      }
      n.double.bonds <-  (rdb.equiv - 3)
      n.fa.carbons <- (c - 9)
      if(n.double.bonds > max.dbl.bnds){
        return(NA)
      }
      return(paste0("MGDG(",n.fa.carbons,":",n.double.bonds,")"))
    }
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
        #some FFA 
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
        }
        if(o == 4){
          return(paste0("phytoCer(",n.fa.carbons,":",n.double.bonds,")"))
        }
        if(o == 5){
          return(paste0("phytoCer(h2,",n.fa.carbons,":",n.double.bonds,")"))
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
          if((c %% 2) == 0){
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
assign_comp_from_mz <- function(mz,elements,isotope,atomic.mass,mins.list,maxs.list,rdb.rule,rdbrange,ionmode,all.elements,error.ppm){
  mz <- as.numeric(mz)
  if(length(mz)>1){
    multiple.mzs <- TRUE
    l<-length(mz)
  }else{
    multiple.mzs <- FALSE
  }
  
  atomic.mass <- as.numeric(atomic.mass)
  mins.list <- as.numeric(mins.list)
  maxs.list <- as.numeric(maxs.list)
  
  if(rdb.rule == "none"){
    v <- c(atomic.mass,1/1837)
    ele.list <- c(mapply(paste0,elements,isotope),"e")
    
    if(ionmode == "neutral"){
      e.bound = 0
    }
    if(ionmode == "neg.ion"){
      e.bound = 1
    }
    if(ionmode == "pos.ion"){
      e.bound = -1
    }

    bounds <- list(lower = list(ind = c(1:length(v)), val = c(mins.list,e.bound)),
                   upper = list(ind = c(1:length(v)), val = c(maxs.list,e.bound)))
    
    rdb.row <- rep(0,length(v))
    rdb.row[ele.list %in% c("H1","D2","T3","Cl35","Cl37")] <- (-0.5)
    rdb.row[ele.list %in% c("C12","C13","C14")] <- 1
    rdb.row[ele.list %in% c("N14","N15","P31")] <- 0.5
    if(error.ppm == "" | is.null(error.ppm) | is.na(error.ppm)){
      error.ppm <- 5
    }
    
    clean_elements <- ifelse(paste0(elements,"[",isotope,"]") %in% all.elements[all.elements$isotopic.composition > 0.9,colnames(all.elements) == "full.element"],elements,paste0(elements,"[",isotope,"]"))
    order <- c("C","C[13]","C[14]","H","D[2]","T[3]","O","O[17]","O[18}","N","N[15]","P","Na","Cl[35]","Cl[37]")
    
    return(
      lapply(mz,
             function(x) {
               if(multiple.mzs){(incProgress(1/l))}
               a<-Rglpk::Rglpk_solve_LP(obj = v, mat = rbind(t(v),rdb.row,rdb.row), dir = c("<=","<=",">="), rhs = c(x,rdbrange[2]-1,rdbrange[1]-1), bounds = bounds,max = TRUE, types = "I")
               b<-Rglpk::Rglpk_solve_LP(obj = v, mat = rbind(t(v),rdb.row,rdb.row), dir = c(">=","<=",">="), rhs = c(x,rdbrange[2]-1,rdbrange[1]-1), bounds = bounds,max = FALSE, types = "I")
               if(((abs(a$optimum-x)/a$optimum)*10^6) <= error.ppm | ((abs(b$optimum-x)/b$optimum)*10^6) <= error.ppm){
                 if(abs(a$optimum-x) < abs(b$optimum-x)){
                   tempdf<-a$solution[!(ele.list %in% c("e","rdbvar"))]
                   clean_elements <- clean_elements[tempdf != 0]
                   tempdf <- tempdf[tempdf != 0]
                   tempdf[tempdf == 1] <- ""
                   return(
                     c(paste(mapply(paste0,
                                    clean_elements[ match(c(order[order%in%clean_elements],clean_elements[!(clean_elements%in%order)]),clean_elements) ],
                                    tempdf[ match(c(order[order%in%clean_elements],clean_elements[!(clean_elements%in%order)]),clean_elements) ]),
                             collapse = " "),
                       a$optimum,
                       ((x-a$optimum)/x)*10^6))
                 }else{
                   tempdf<-b$solution[!(ele.list %in% c("e","rdbvar"))]
                   clean_elements <- clean_elements[tempdf != 0]
                   tempdf <- tempdf[tempdf != 0]
                   tempdf[tempdf == 1] <- ""
                   return(
                     c(paste(mapply(paste0,
                                    clean_elements[ match(c(order[order%in%clean_elements],clean_elements[!(clean_elements%in%order)]),clean_elements) ],
                                    tempdf[ match(c(order[order%in%clean_elements],clean_elements[!(clean_elements%in%order)]),clean_elements) ]),
                             collapse = " "),
                       b$optimum,
                       ((b$optimum-x)/x)*10^6))
                 }
                 }else{
                   return(c(rep("",3)))
                 }
               }
             )
      ) 
    
    
  }else{
    if(rdb.rule == "nonint"){
      rhs.var = 0.5
    }else{
      #rdb rule = int
      rhs.var = 0
    }
    
    v <- c(atomic.mass,1/1837,0)
    ele.list <- c(mapply(paste0,elements,isotope),"e","rdbvar")
    
    if(ionmode == "neutral"){
      e.bound = 0
    }
    if(ionmode == "neg.ion"){
      e.bound = 1
    }
    if(ionmode == "pos.ion"){
      e.bound = -1
    }
    
    min.r = floor((rhs.var-rdbrange[2])/2)-1
    max.r = ceiling((rhs.var-rdbrange[1])/2)+1

    bounds <- list(lower = list(ind = c(1:length(v)), val = c(mins.list,e.bound,min.r)),
                   upper = list(ind = c(1:length(v)), val = c(maxs.list,e.bound,max.r)))
    
    rdb.row <- rep(0,length(v))
    rdb.row[ele.list %in% c("H1","D2","T3","Cl35","Cl37")] <- (-0.5)
    rdb.row[ele.list %in% c("C12","C13","C14")] <- 1
    rdb.row[ele.list %in% c("N14","N15","P31")] <- 0.5
    
    rdb.row2 <- rdb.row
    rdb.row2[length(rdb.row2)] <- 1
    
    if(error.ppm == "" | is.null(error.ppm) | is.na(error.ppm)){
      error.ppm <- 5
    }
    
    clean_elements <- ifelse(paste0(elements,"[",isotope,"]") %in% all.elements[all.elements$isotopic.composition > 0.9,colnames(all.elements) == "full.element"],elements,paste0(elements,"[",isotope,"]"))
    order <- c("C","C[13]","C[14]","H","D[2]","T[3]","O","O[17]","O[18}","N","N[15]","P","Na","Cl[35]","Cl[37]")
    
    
    return(lapply(mz,function(x){
      if(multiple.mzs){(incProgress(1/l))}
    a<-Rglpk::Rglpk_solve_LP(obj = v, mat = rbind(t(v),rdb.row,rdb.row,rdb.row2), dir = c("<=","<=",">=","=="), rhs = c(x,rdbrange[2]-1,rdbrange[1]-1,rhs.var), bounds,max = TRUE, types = "I")
    b<-Rglpk::Rglpk_solve_LP(obj = v, mat = rbind(t(v),rdb.row,rdb.row,rdb.row2), dir = c(">=","<=",">=","=="), rhs = c(x,rdbrange[2]-1,rdbrange[1]-1,rhs.var), bounds,max = FALSE, types = "I")

    if(((abs(a$optimum-x)/a$optimum)*10^6) <= error.ppm | ((abs(b$optimum-x)/b$optimum)*10^6) <= error.ppm){
      if(abs(a$optimum-x) < abs(b$optimum-x)){
        tempdf<-a$solution[!(ele.list %in% c("e","rdbvar"))]
        clean_elements <- clean_elements[tempdf != 0]
        tempdf <- tempdf[tempdf != 0]
        tempdf[tempdf == 1] <- ""
        return(c(paste(mapply(paste0,
                              clean_elements[ match(c(order[order%in%clean_elements],clean_elements[!(clean_elements%in%order)]),clean_elements) ],
                              tempdf[ match(c(order[order%in%clean_elements],clean_elements[!(clean_elements%in%order)]),clean_elements) ]),
                       collapse = " "),
                 a$optimum,
                 ((x-a$optimum)/x)*10^6))
      }else{
        tempdf<-b$solution[!(ele.list %in% c("e","rdbvar"))]
        clean_elements <- clean_elements[tempdf != 0]
        tempdf <- tempdf[tempdf != 0]
        tempdf[tempdf == 1] <- ""
        return(c(paste(mapply(paste0,
                              clean_elements[ match(c(order[order%in%clean_elements],clean_elements[!(clean_elements%in%order)]),clean_elements) ],
                              tempdf[ match(c(order[order%in%clean_elements],clean_elements[!(clean_elements%in%order)]),clean_elements) ]),
                       collapse = " "),
                 b$optimum,
                 ((b$optimum-x)/x)*10^6))
      }
    }else{
      return(c(rep("",3)))
    }}))
    
  }
  
}
calculate_logical_element_minimums <- function(loi,lipidclasslist,elements,isotopes,current.mins){
  ############
  min.c.12 <- c(33, #dag
             43, #tag
             15, #pa
             15, #pg
             25, #pi
             40, #cl
             15, #sm
             15, #cer
             15, #pc
             15, #pe
             15, #ps
             13, #ffa
             30 #gdg
             )

  min.h.1 <- c(30 #for all
               )

  min.o.16 <- c(4, #dag
                5, #tag
                7, #pa
                9, #pg
                12, #pi
                16, #cl
                5, #sm
                3, #cer
                7, #pc
                7, #pe
                9, #ps
                2, #ffa
                10 #gdg
  )

  min.n.14 <- c(0, #dag
                0, #tag
                0, #pa
                0, #pg
                0, #pi
                0, #cl
                2, #sm
                1, #cer
                1, #pc
                1, #pe
                1, #ps
                0, #ffa
                0 #gdg
  )

  min.p.31 <- c(0, #dag
                0, #tag
                1, #pa
                1, #pg
                0, #pi
                2, #cl
                1, #sm
                0, #cer
                1, #pc
                1, #pe
                1, #ps
                0, #ffa
                0 #gdg
  )

  ##############
  mins <- current.mins 
  
  mins[paste0(elements,isotopes) == "C12"] <- min(min.c.12[lipidclasslist %in% loi]) 
  mins[paste0(elements,isotopes) == "H1"] <- min.h.1[1] 
  mins[paste0(elements,isotopes) == "O16"] <- min(min.o.16[lipidclasslist %in% loi]) 
  mins[paste0(elements,isotopes) == "N14"] <- min(min.n.14[lipidclasslist %in% loi])
  mins[paste0(elements,isotopes) == "P31"] <- min(min.p.31[lipidclasslist %in% loi]) 
  
  return(mins)
  
}
calculate_logical_element_maximums <- function(loi,lipidclasslist,elements,isotopes,current.maxs){
  ############
  max.c.12 <- c(65, #dag
                75, #tag
                75, #pa
                75, #pg
                75, #pi
                90, #cl
                75, #sm
                75, #cer
                75, #pc
                75, #pe
                75, #ps
                65, #ffa
                70 #gdg
  )

  max.h.1 <- c(120, #dag
               130, #tag
               120, #pa
               120, #pg
               120, #pi
               160, #cl
               120, #sm
               120, #cer
               120, #pc
               120, #pe
               120, #ps
               120, #ffa
               160 #gdg
  )

  max.o.16 <- c(6, #dag
                7, #tag
                8, #pa
                10, #pg
                13, #pi
                17, #cl
                7, #sm
                4, #cer
                8, #pc
                8, #pe
                10, #ps
                4, #ffa
                15 #gdg
  )

  max.n.14 <- c(1, #dag
                1, #tag
                0, #pa
                0, #pg
                0, #pi
                0, #cl
                2, #sm
                1, #cer
                1, #pc
                1, #pe
                1, #ps
                0, #ffa
                0 #gdg
  )

  max.p.31 <- c(0, #dag
                0, #tag
                1, #pa
                1, #pg
                1, #pi
                2, #cl
                1, #sm
                0, #cer
                1, #pc
                1, #pe
                1, #ps
                0, #ffa
                0 #gdg
  )
  ##############
  maxs <- current.maxs 
  
  maxs[paste0(elements,isotopes) == "C12"] <- max(max.c.12[lipidclasslist %in% loi]) 
  maxs[paste0(elements,isotopes) == "H1"] <- max(max.h.1[lipidclasslist %in% loi]) 
  maxs[paste0(elements,isotopes) == "O16"] <- max(max.o.16[lipidclasslist %in% loi]) 
  maxs[paste0(elements,isotopes) == "N14"] <- max(max.n.14[lipidclasslist %in% loi])
  maxs[paste0(elements,isotopes) == "P31"] <- max(max.p.31[lipidclasslist %in% loi]) 
  
  return(maxs)
  
}

# make lists
lipidclasslist <- c("DAG" = "dag","TAG" = "tag","PA" = "pa","PG" = "pg","PI" = "pi","CL" = "cl","SM" = "sm","Ceramide" = "cer","PC" = "pc","PE" = "pe","PS" = "ps","FFA" = "ffa","GDG" = "gdg")
negionclasses <- c("pe","cer","cl","pi","pg","pa","ps","ffa","gdg")
posionclasses <- c("pc","sm","tag","dag","gdg")
relevant.col.names <- c("mz","intensity","rel.intensity","theo.mass","delta","rdb.equiv","composition")

pls <- read.csv('pl_database.csv')
strainslist <- c("n/a",unique(pls$strain))
strainslist <- strainslist[strainslist != ""]
names(strainslist) <- strainslist

all.elements = read.csv('atomics_wts_and_isotopic_compositions_for_all_elements.csv')
all.elements$full.element = paste0(all.elements$element,"[",all.elements$isotope,"]")
this_table = all.elements[paste0(all.elements$element,all.elements$isotope) %in% c("H1","C12","N14","P31","O16"),
                          colnames(all.elements) %in% c("element","isotope","atomic.mass")]
this_table$minimum <- rep(0,nrow(this_table))
this_table$maximum <- rep(150,nrow(this_table))
rownames(this_table) <- NULL
this_table2 <- this_table
choices = paste0(all.elements$element,"[",all.elements$isotope,"]")
names(choices) = all.elements$full.element

ui <- fluidPage(
  titlePanel("SLAT: Shotgun Lipidomic Assignment Tool"),
  tabsetPanel(
    id = "tabsid",
    tabPanel(value = 0,
             titlePanel("Analyze a single peak"),
             fluidRow(
               column(width = 3,
                 h1(strong("1. Enter general analysis parameters:"), style = "font-size:22.5px;"),
                 selectInput("domain","Sample source",
                             c("Eukaryota" = "euk",
                               "Bacteria" = "bact"),
                             selected = "euk"),
                 conditionalPanel(
                   condition = "input.tabsid == 0  && input.domain == 'bact'",
                   selectInput("strain","Bacterial strain:",
                               strainslist,
                               selected = "n/a",
                               multiple = F)
                 ),
                 selectInput("ionmode", "What mode were samples run in?",
                             c("negative ion mode" = "neg.ion",
                               "positive ion mode" = "pos.ion",
                               "N/A (search for structure as a neutral species)" = "neutral")),
                 conditionalPanel(
                   condition = "input.tabsid == 0  && input.ionmode == 'pos.ion'",
                   checkboxGroupInput("posionadduct", "Ions/adducts to look for:",
                                      c("[M+H]+" = "m.plus.h",
                                        "[M+NH4]+" = "m.plus.ammonia",
                                        "[M+Na]+" = "m.plus.sodium")),
                   selectInput("posionloi","Lipid classes to search for:",
                               lipidclasslist,
                               selected = posionclasses,
                               multiple = T)
                 ),
                 conditionalPanel(
                   condition = "input.tabsid == 0  && input.ionmode == 'neg.ion'",
                   checkboxGroupInput("negionadduct", "Ions/adducts to look for:",
                                      c("[M-H]-" = "m.minus.h",
                                        "[M+Cl]-" = "m.plus.chloride",
                                        "[M-2H]2-" = "m.minus.2h")),
                   selectInput("negionloi","Lipid classes to search for:",
                               lipidclasslist,
                               selected = negionclasses,
                               multiple = T)
                 ),
                 conditionalPanel(
                   condition = "input.tabsid == 0  && input.ionmode == 'neutral'",
                   selectInput("neutralloi","Lipid classes to search for:",
                               lipidclasslist,
                               selected = c(negionclasses,posionclasses),
                               multiple = T)
                 ),
                 sliderInput("maxdblbnds",
                             "Maximum number of acyl/alkyl double bonds allowed",
                             min = 0,max = 16,value = 12,step = 1,round = T,animate = F, ticks = F),
                 h1(strong("2. Enter a composition or m/z value:"), style = "font-size:22.5px;"),
                 textInput("comp", label = NULL, value = "", width = NULL, placeholder = "e.g. C36 H73 O8 N P or 634.4452")
               ), #end first column 
               column(width = 6,
                      conditionalPanel(
                        condition = "input.tabsid == 0  && !isNaN(input.comp) && input.comp != '' ",
                        h1(strong("3. Set parameters for assignment of composition from m/z:"), style = "font-size:22.5px;"),
                        splitLayout(radioButtons("rdbrule","Force RDB to equal:", 
                                                 choices = c("non-integer" = "nonint",
                                                             "integer" = "int",
                                                             "N/A (no rule)" = "none"),
                                                 selected = "none",inline = T),
                                    sliderInput("rdbrange",
                                                "RDB range:",
                                                min = -50, 
                                                max = 100, 
                                                value = c(-1,100),step = 0.5)),
                        selectInput("element",
                                                label = "Element:",
                                                choices = choices,
                                                selected = "C[12]"),
                        textInput("atomic.mass",
                                              "Atomic mass:",
                                              value = 12.0000000),
                        splitLayout(
                          textInput(inputId="ele.min", label="Element min:", value = 20),
                          textInput(inputId="ele.max", label="Element max:", value = 120)
                        ),
                        numericInput("delta.error",
                                     "Mass tolerance (ppm):",
                                      value = 3,
                                      min = 0.5, max = 100, step = 0.5),
                        actionButton("add_btn", "Add/Update element"),
                        actionButton("delete_btn", "Delete element"),
                        DTOutput("shiny_table",width = "20%")
                      ),
                      conditionalPanel(
                        condition = "input.tabsid == 0  && (isNaN(input.comp) || input.comp == '')",
                        htmlOutput("analyzed_data1")
                      )
               ),#end second column
               column(width = 3,
                        conditionalPanel(condition = "input.tabsid == 0  && !isNaN(input.comp) && input.comp != '' ",
                                         htmlOutput("analyzed_data2")
                                         )
                      )
               )
    ),
    
    tabPanel(value = 1,
             titlePanel("Paste all data for a sample"),
             fluidRow(
               column(
                 width = 3,
                 h1(strong("1. Enter general parameters:"), style = "font-size:22.5px;"),
                 selectInput("bulkdomain","Sample source",
                             c("Eukaryota" = "euk",
                               "Bacteria" = "bact"),
                             selected = "euk"),
                 conditionalPanel(
                   condition = "input.tabsid == 1  && input.bulkdomain == 'bact'",
                   selectInput("bulkstrain","Bacterial strain:",
                               strainslist,
                               selected = "n/a",
                               multiple = F)
                 ),
                 selectInput("bulkionmode", "What mode were samples run in?",
                             c("negative ion mode" = "neg.ion",
                               "positive ion mode" = "pos.ion",
                               "N/A (search for structure as a neutral species)" = "neutral")) ,
                 conditionalPanel(
                   condition = "input.tabsid == 1  && input.bulkionmode == 'pos.ion'",
                   checkboxGroupInput("bulkposionadduct", "Ions/adducts to look for:",
                                      c("[M+H]+" = "m.plus.h",
                                        "[M+NH4]+" = "m.plus.ammonia",
                                        "[M+Na]+" = "m.plus.sodium")),
                   selectInput("bulkposionloi","Lipid classes to search for:",
                               lipidclasslist,
                               selected = posionclasses,
                               multiple = T)
                 ),
                 conditionalPanel(
                   condition = "input.tabsid == 1  && input.bulkionmode == 'neg.ion'",
                   checkboxGroupInput("bulknegionadduct", "Ions/adducts to look for:",
                                      c("[M-H]-" = "m.minus.h",
                                        "[M+Cl]-" = "m.plus.chloride",
                                        "[M-2H]2-" = "m.minus.2h")),
                   selectInput("bulknegionloi","Lipid classes to search for:",
                               lipidclasslist,
                               selected = negionclasses,
                               multiple = T)
                 ),
                 conditionalPanel(
                   condition = "input.tabsid == 1  && input.bulkionmode == 'neutral'",
                   selectInput("bulkneutralloi","Lipid classes to search for:",
                               lipidclasslist,
                               selected = c(negionclasses,posionclasses),
                               multiple = T)
                 ),
                 radioButtons("returnAll", "Would you like to...", 
                              choiceNames = c("Return all peaks",
                                              "Only return peaks with a species assigned"),
                              choiceValues = c("returnallpeaks","returnassignedpeaks"),
                              inline = F, selected = ("returnallpeaks")),
                 sliderInput("bulkmaxdblbnds",
                             "Max number of acyl/alkyl double bonds allowed",
                             min = 0,max = 16,value = 12,step = 1,round = T,animate = F, ticks = F),
                 h1(strong("2. Select starting point:"), style = "font-size:22.5px;"),
                 radioButtons("assignmentmode", "Assign starting from m/z values or pre-determined compositions?", 
                              choiceNames = c("From m/z values",
                                              "From compositions"),
                              choiceValues = c("frommz","fromcomp"),
                              inline = T, selected = character(0)) 
                 ),
               column(width = 5,
                      conditionalPanel(
                        condition = "input.tabsid == 1  && input.assignmentmode == 'fromcomp'",
                        h1(strong("3. Enter a few more parameters:"), style = "font-size:22.5px;"),
                        splitLayout(selectInput("bulkisotope", "Isotope being used",
                                                c("None" = "none",
                                                  "[13]C" = "carbon_isotope",
                                                  "[2]H" = "hydrogen_isotope")),
                        numericInput("deltacutoff", "abs(Delta) cutoff:", 10, min = 0.1, max = 1000,value = 1000)),
                        h1(strong("4. Input data:"), style = "font-size:22.5px;"),
                        import_copypaste_ui2("myid",
                                             title = "")
                      ),
                      tags$head(tags$style(".shiny-notification {position: fixed; top: 50% ;left: 35%; color: red; font-weight: bold; font-size: 30px;} .progress-bar {background-color: red;}")),
                      conditionalPanel(condition = "input.tabsid == 1  && input.assignmentmode == 'frommz'",
                                       h1(strong("3. Set parameters for assignment of composition from m/z:"), style = "font-size:22.5px;"),
                                       splitLayout(radioButtons("bulkrdbrule","Force RDB to equal:", 
                                                                choices = c("non-integer" = "nonint",
                                                                            "integer" = "int",
                                                                            "N/A (no rule)" = "none"),
                                                                selected = "none",inline = T),
                                                   sliderInput("bulkrdbrange",
                                                               "RDB range:",
                                                               min = -50, 
                                                               max = 100, 
                                                               value = c(-1,100),step = 0.5,ticks = F)),
                                       splitLayout(
                                         selectInput("bulkelement",
                                                     label = "Element:",
                                                     choices = choices,
                                                     selected = "C[12]"),
                                         textInput("bulkatomic.mass",
                                                     "Atomic mass:",
                                                     value = 12.0000000),
                                         tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible;}"))),
                                         cellWidths = c("50%", "50%")),
                                       splitLayout(
                                         textInput(inputId="bulk.ele.min", label="Element min:", value = 20),
                                         textInput(inputId="bulk.ele.max", label="Element max:", value = 120)
                                       ),
                                       actionButton("bulk.add_btn", "Add/Update element"),
                                       actionButton("bulk.delete_btn", "Delete element"),
                                       DTOutput("bulk_shiny_table",width = "20%"),
                                       conditionalPanel(condition = "input.tabsid == 1 && input.transferdata == 'y' ",
                                                        h1(strong("4. Recalibrated data has been automatically imported from the recalibration tab"), style = "font-size:22.5px;")),
                                       conditionalPanel(condition = "input.tabsid == 1 && input.transferdata !== 'y' ",
                                                        h1(strong("4. Input data:"), style = "font-size:22.5px;"),
                                                        import_copypaste_ui3("myid2",
                                                                             title = "")
                                                        )
                                       )
                      ),
               column(
                 width = 4,
                 conditionalPanel(condition = "input.tabsid == 1 && input.assignmentmode == 'fromcomp'",
                                  h1(strong("5. Preview results:"), style = "font-size:22.5px;"),
                                  tags$b("Import status:"),
                                  htmlOutput(outputId = "status"),
                                  tags$b("Sample name:"),
                                  verbatimTextOutput(outputId = "name"),
                                  tags$b("Verify that data was input & read correctly (just make sure columns align):"),
                                  verbatimTextOutput(outputId = "data"),
                                  h1(strong("6. Download results:"), style = "font-size:22.5px;"),
                                  downloadButton("download", "Download results")),
                 conditionalPanel(condition = "input.tabsid == 1 && input.assignmentmode == 'frommz'",
                                  h1(strong("5. Set tolerance and analyze data:"), style = "font-size:22.5px;"),
                                  splitLayout(
                                    numericInput("bulk.delta.error",
                                                 "Mass tolerance (ppm):",
                                                 value = 3,
                                                 min = 0.5, max = 100, step = 0.5),
                                    actionButton("bulk.analyze.button", "Analyze data", icon("redo"),style="color: #fff; background-color: #FF0000; border-color: #2e6da4; margin-top:25px",width = '100%')
                                  ),
                                  tags$b("Import status:"),
                                  htmlOutput(outputId = "status2"),
                                  tags$b("Sample name:"),
                                  verbatimTextOutput(outputId = "name2"),
                                  tags$b("Verify that data was input & read correctly:"),
                                  verbatimTextOutput(outputId = "data2"),
                                  h1(strong("6. Download results:"), style = "font-size:22.5px;"),
                                  downloadButton("download2", "Download results")
                                  )
                 )
               ),
             ),
    tabPanel(value = 2,
             titlePanel("Recalibrate raw data"),
             fluidRow(
               column(
                 width = 4,
                 h1(strong("1. Input data:"), style = "font-size:22.5px;"),
                 import_copypaste_ui3("myid3",
                                      title = ""),
                 tags$b("Import status:"),
                 htmlOutput(outputId = "status3"),
                 tags$b("Sample name:"),
                 verbatimTextOutput(outputId = "name3")
               ),
               column(width = 8,
                      h1(strong("2. Click and drag to zoom in on an m/z range, click a point to select it, or double click to reset"), style = "font-size:22.5px;"),
                      plotlyOutput("plot1"),
                      h1(strong("3. Input exact theoretical mass of selected points:"), style = "font-size:22.5px;"),
                      splitLayout(DTOutput("alignment_shiny_table",width = "50%"),
                                  actionButton("alignment_delete_btn", "Delete selected point",icon("trash"),style="color: #FFFFFF; background-color: #FF0000; border-color: #2e6da4; margin-top:25px")),
                      h1(strong("4. Copy or download results:"), style = "font-size:22.5px;"),
                      rclipboardSetup(),
                      splitLayout(radioButtons("transferdata",
                                               "Transfer recalibrated data to second tab?",
                                               choices = c("Yes" = "y",
                                                           "No" = "n"),
                                               selected = character(0),
                                               inline = T),
                                  uiOutput("clip"),
                                  downloadButton("download3", "Download aligned data")
                                  )
                      )
             )),
    )
  
)

server <- function(input, output, session){
  observeEvent(input$negionadduct,{
    if("m.plus.chloride" %in% input$negionadduct){
      f.e <- "Cl[35]"
      t <- data.frame(element = all.elements$element[match(f.e,all.elements$full.element)],
                      isotope = all.elements$isotope[match(f.e,all.elements$full.element)],
                      atomic.mass = all.elements$atomic.mass[match(f.e,all.elements$full.element)],
                      minimum = 0,
                      maximum = 1)
      if( !(paste0(t$element,t$isotope) %in% paste0(this_table()$element,this_table()$isotope)) ){
        this_table(rbind(t, this_table()))
      }
    }else{
      this_table(this_table()[this_table()$element != "Cl",])
    }
  })
  observeEvent(input$bulknegionadduct,{
    if("m.plus.chloride" %in% input$bulknegionadduct){
      f.e <- "Cl[35]"
      t <- data.frame(element = all.elements$element[match(f.e,all.elements$full.element)],
                      isotope = all.elements$isotope[match(f.e,all.elements$full.element)],
                      atomic.mass = all.elements$atomic.mass[match(f.e,all.elements$full.element)],
                      minimum = 0,
                      maximum = 1)
      if( !(paste0(t$element,t$isotope) %in% paste0(this_table2()$element,this_table2()$isotope)) ){
        this_table2(rbind(t, this_table2()))
      }
    }else{
      this_table2(this_table2()[this_table2()$element != "Cl",])
    }
  })
  observeEvent(input$posionadduct,{
    if("m.plus.sodium" %in% input$posionadduct){
      f.e <- "Na[23]"
      t <- data.frame(element = all.elements$element[match(f.e,all.elements$full.element)],
                      isotope = all.elements$isotope[match(f.e,all.elements$full.element)],
                      atomic.mass = all.elements$atomic.mass[match(f.e,all.elements$full.element)],
                      minimum = 0,
                      maximum = 1)
      if( !(paste0(t$element,t$isotope) %in% paste0(this_table()$element,this_table()$isotope)) ){
        this_table(rbind(t, this_table()))
      }
    }else{
      this_table(this_table()[this_table()$element != "Na",])
    }
  })
  observeEvent(input$bulkposionadduct,{
    if("m.plus.sodium" %in% input$bulkposionadduct){
      f.e <- "Na[23]"
      t <- data.frame(element = all.elements$element[match(f.e,all.elements$full.element)],
                      isotope = all.elements$isotope[match(f.e,all.elements$full.element)],
                      atomic.mass = all.elements$atomic.mass[match(f.e,all.elements$full.element)],
                      minimum = 0,
                      maximum = 1)
      if( !(paste0(t$element,t$isotope) %in% paste0(this_table2()$element,this_table2()$isotope)) ){
        this_table2(rbind(t, this_table2()))
      }
    }else{
      this_table2(this_table2()[this_table2()$element != "Na",])
    }
  })
  #reset strain to n/a everytime domain is changed
  observeEvent(input$domain,
               updateSelectInput(session = session,
                                 "strain",
                                 "Bacterial strain:",
                                 choices = strainslist,
                                 selected = "n/a"
               )
  )
  observeEvent(input$bulkdomain,
               updateSelectInput(session = session,
                                 "bulkstrain",
                                 "Bacterial strain:",
                                 choices = strainslist,
                                 selected = "n/a"
               )
  )
  
  #when user selects a new element from dropdown, autofill its atomic mass in the atomic mass section
  observeEvent(input$element,
               updateTextInput(session = session,
                               "atomic.mass",
                               "Atomic mass:",
                               value = all.elements$atomic.mass[match(input$element,all.elements$full.element)])
  )
  observeEvent(input$bulkelement,
               updateTextInput(session = session,
                               "bulkatomic.mass",
                               "Atomic mass:",
                               value = all.elements$atomic.mass[match(input$bulkelement,all.elements$full.element)])
  )
  
  #make the datatable a reactiveVal
  # this_table is for single datapoint tab, this_table2 is for bulk processing tab
  this_table <- reactiveVal(this_table)
  this_table2 <- reactiveVal(this_table2)
  

  #add element to table when user presses button to do so
  observeEvent(input$add_btn, {
    t <- data.frame(element = all.elements$element[match(input$element,all.elements$full.element)],
                    isotope = all.elements$isotope[match(input$element,all.elements$full.element)],
                    atomic.mass = input$atomic.mass,
                    minimum = ifelse(is.na(input$ele.min) | is.null(input$ele.min) | input$ele.min == "",0,input$ele.min),
                    maximum = ifelse(is.na(input$ele.max) | is.null(input$ele.max) | input$ele.max == "",120,input$ele.max))
    
    if(paste0(t$element,t$isotope) %in% paste0(this_table()$element,this_table()$isotope)){
      c<-this_table()
      c[paste0(c$element,c$isotope) == paste0(t$element,t$isotope),]<-t
      this_table(c)
    }else{
      this_table(rbind(t, this_table()))
    }
  })
  observeEvent(input$bulk.add_btn, {
    t <- data.frame(element = all.elements$element[match(input$bulkelement,all.elements$full.element)],
                    isotope = all.elements$isotope[match(input$bulkelement,all.elements$full.element)],
                    atomic.mass = input$bulkatomic.mass,
                    minimum = ifelse(is.na(input$bulk.ele.min) | is.null(input$bulk.ele.min) | input$bulk.ele.min == "",0,input$bulk.ele.min),
                    maximum = ifelse(is.na(input$bulk.ele.max) | is.null(input$bulk.ele.max) | input$bulk.ele.max == "",120,input$bulk.ele.max))
    
    if(paste0(t$element,t$isotope) %in% paste0(this_table2()$element,this_table2()$isotope)){
      c<-this_table2()
      c[paste0(c$element,c$isotope) == paste0(t$element,t$isotope),]<-t
      this_table2(c)
    }else{
      this_table2(rbind(t, this_table2()))
    }
  })
  
  #display row info when user clicks on a row
  observeEvent(input$shiny_table_rows_selected,
               updateSelectInput(session = session,
                                 "element",
                                 choices = choices,
                                 selected = paste0(this_table()$element[input$shiny_table_rows_selected],"[",this_table()$isotope[input$shiny_table_rows_selected],"]")
                                 )
               )
  observeEvent(input$shiny_table_rows_selected,
               updateTextInput(session = session,
                               "ele.min",
                               label="Element min",
                               value = this_table()$minimum[input$shiny_table_rows_selected])
               )
  observeEvent(input$shiny_table_rows_selected,
               updateTextInput(session = session,
                               "ele.max",
                               label="Element max",
                               value = this_table()$maximum[input$shiny_table_rows_selected])
               )
  observeEvent(input$bulk_shiny_table_rows_selected,
               updateSelectInput(session = session,
                                 "bulkelement",
                                 choices = choices,
                                 selected = paste0(this_table2()$element[input$bulk_shiny_table_rows_selected],"[",this_table2()$isotope[input$bulk_shiny_table_rows_selected],"]")
               )
  )
  observeEvent(input$bulk_shiny_table_rows_selected,
               updateTextInput(session = session,
                               "bulk.ele.min",
                               label="Element min",
                               value = this_table2()$minimum[input$bulk_shiny_table_rows_selected])
  )
  observeEvent(input$bulk_shiny_table_rows_selected,
               updateTextInput(session = session,
                               "bulk.ele.max",
                               label="Element max",
                               value = this_table2()$maximum[input$bulk_shiny_table_rows_selected])
  )
  
  #reset elemental ranges to logical max/mins everytime a lipid is added or removed
  observeEvent(input$negionloi,{
    t<-this_table()
    t$minimum <- calculate_logical_element_minimums(input$negionloi,lipidclasslist,t$element,t$isotope,t$minimum)
    t$maximum <- calculate_logical_element_maximums(input$negionloi,lipidclasslist,t$element,t$isotope,t$maximum)
    this_table(t)
  })
  observeEvent(input$bulknegionloi,{
    t<-this_table2()
    t$minimum <- calculate_logical_element_minimums(input$bulknegionloi,lipidclasslist,t$element,t$isotope,t$minimum)
    t$maximum <- calculate_logical_element_maximums(input$bulknegionloi,lipidclasslist,t$element,t$isotope,t$maximum)
    this_table2(t)
  })
  observeEvent(input$posionloi,{
    t<-this_table()
    t$minimum <- calculate_logical_element_minimums(input$posionloi,lipidclasslist,t$element,t$isotope,t$minimum)
    t$maximum <- calculate_logical_element_maximums(input$posionloi,lipidclasslist,t$element,t$isotope,t$maximum)
    this_table(t)
  })
  observeEvent(input$bulkposionloi,{
    t<-this_table2()
    t$minimum <- calculate_logical_element_minimums(input$bulkposionloi,lipidclasslist,t$element,t$isotope,t$minimum)
    t$maximum <- calculate_logical_element_maximums(input$bulkposionloi,lipidclasslist,t$element,t$isotope,t$maximum)
    this_table2(t)
  })
  observeEvent(input$neutralloi,{
    t<-this_table()
    t$minimum <- calculate_logical_element_minimums(input$neutralloi,lipidclasslist,t$element,t$isotope,t$minimum)
    t$maximum <- calculate_logical_element_maximums(input$neutralloi,lipidclasslist,t$element,t$isotope,t$maximum)
    this_table(t)
  })
  observeEvent(input$bulkneutralloi,{
    t<-this_table2()
    t$minimum <- calculate_logical_element_minimums(input$bulkneutralloi,lipidclasslist,t$element,t$isotope,t$minimum)
    t$maximum <- calculate_logical_element_maximums(input$bulkneutralloi,lipidclasslist,t$element,t$isotope,t$maximum)
    this_table2(t)
  })

  #delete a row (= element) when user clicks button to do so
  observeEvent(input$delete_btn, {
    t = this_table()
    if (!is.null(input$shiny_table_rows_selected)) {
      t <- t[-as.numeric(input$shiny_table_rows_selected),]
    }
    this_table(t)
  })
  observeEvent(input$bulk.delete_btn, {
    t = this_table2()
    if (!is.null(input$bulk_shiny_table_rows_selected)) {
      t <- t[-as.numeric(input$bulk_shiny_table_rows_selected),]
    }
    this_table2(t)
  })
  
  #render datatable
  output$shiny_table <- renderDT({
    datatable(this_table(), selection = 'single', options = list(dom = 't'))
  })
  output$bulk_shiny_table <- renderDT({
    datatable(this_table2(), selection = 'single', options = list(dom = 't'))
  })
  
  #every time user changes ion mode, reset ion and lipid selections for all other ion modes than the one they changed selection to
  observeEvent(input$ionmode,
               {if(input$ionmode == "pos.ion"){
                 #reset negative and neutral
                 updateCheckboxGroupInput(session = session,
                                          "negionadduct",
                                          "Ions/adducts to look for:",
                                          c("[M-H]-" = "m.minus.h",
                                            "[M+Cl]-" = "m.plus.chloride",
                                            "[M-2H]2-" = "m.minus.2h"),
                                          selected=character(0))
                 updateSelectInput(session = session, "negionloi",
                                   selected = negionclasses)
                 updateSelectInput(session = session, "neutralloi",
                                   selected = c(negionclasses,posionclasses))
                 t<-this_table()
                 t$minimum <- calculate_logical_element_minimums(input$posionloi,lipidclasslist,t$element,t$isotope,t$minimum)
                 t$maximum <- calculate_logical_element_maximums(input$posionloi,lipidclasslist,t$element,t$isotope,t$maximum)
                 this_table(t)
                 
               }
                 if(input$ionmode == "neg.ion"){
                   updateCheckboxGroupInput(session = session,
                                            "posionadduct",
                                            "Ions/adducts to look for:",
                                            c("[M+H]+" = "m.plus.h",
                                              "[M+NH4]+" = "m.plus.ammonia",
                                              "[M+Na]+" = "m.plus.sodium"),
                                            selected=character(0))
                   updateSelectInput(session = session, "posionloi",
                                     selected = posionclasses)
                   updateSelectInput(session = session, "neutralloi",
                                     selected = c(negionclasses,posionclasses))
                   t<-this_table()
                   t$minimum <- calculate_logical_element_minimums(input$negionloi,lipidclasslist,t$element,t$isotope,t$minimum)
                   t$maximum <- calculate_logical_element_maximums(input$negionloi,lipidclasslist,t$element,t$isotope,t$maximum)
                   this_table(t)
                   
                 }
                 if(input$ionmode == "neutral"){
                   updateCheckboxGroupInput(session = session,
                                            "negionadduct",
                                            "Ions/adducts to look for:",
                                            c("[M-H]-" = "m.minus.h",
                                              "[M+Cl]-" = "m.plus.chloride",
                                              "[M-2H]2-" = "m.minus.2h"),
                                            selected=character(0))
                   updateSelectInput(session = session, "negionloi",
                                     selected = negionclasses)
                   updateCheckboxGroupInput(session = session,
                                            "posionadduct",
                                            "Ions/adducts to look for:",
                                            c("[M+H]+" = "m.plus.h",
                                              "[M+NH4]+" = "m.plus.ammonia",
                                              "[M+Na]+" = "m.plus.sodium"),
                                            selected=character(0))
                   updateSelectInput(session = session, "posionloi",
                                     selected = posionclasses)
                   t<-this_table()
                   t$minimum <- calculate_logical_element_minimums(input$neutralloi,lipidclasslist,t$element,t$isotope,t$minimum)
                   t$maximum <- calculate_logical_element_maximums(input$neutralloi,lipidclasslist,t$element,t$isotope,t$maximum)
                   this_table(t)
                 }
                 
               },
               ignoreInit = T)
  observeEvent(input$bulkionmode,
               {if(input$bulkionmode == "pos.ion"){
                 updateCheckboxGroupInput(session = session,
                                          "bulknegionadduct",
                                          "Ions/adducts to look for:",
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
                                            "Ions/adducts to look for:",
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
                                            "Ions/adducts to look for:",
                                            c("[M-H]-" = "m.minus.h",
                                              "[M+Cl]-" = "m.plus.chloride",
                                              "[M-2H]2-" = "m.minus.2h"),
                                            selected=character(0))
                   updateSelectInput(session = session, "bulknegionloi",
                                     selected = negionclasses)
                   updateCheckboxGroupInput(session = session,
                                            "bulkposionadduct",
                                            "Ions/adducts to look for:",
                                            c("[M+H]+" = "m.plus.h",
                                              "[M+NH4]+" = "m.plus.ammonia",
                                              "[M+Na]+" = "m.plus.sodium"),
                                            selected=character(0))
                   updateSelectInput(session = session, "bulkposionloi",
                                     selected = posionclasses)
                 }
               },
               ignoreInit = T)
  observeEvent(tab3data(),{
    if(!is.null(imported3$status())){
      if(imported3$status() == "success"){
        updateRadioButtons(session = session,
                           inputId = "assignmentmode",
                           
                           choiceNames = c("From m/z values",
                                           "From compositions"),
                           choiceValues = c("frommz","fromcomp"),
                           inline = T, 
                           selected = "frommz"
                           )
      }
    }
    
  })
  
  #analyze data for single ion mode
  analyzed_data <- renderUI({
    if(!is.na(suppressWarnings(as.numeric(input$comp)))){
      #input is numeric
      #assign composition to m/z:
      comp.object <- assign_comp_from_mz(mz = input$comp,
                                         elements = this_table()$element,
                                         isotope = this_table()$isotope,
                                         atomic.mass = this_table()$atomic.mass,
                                         mins.list = this_table()$minimum,
                                         maxs.list = this_table()$maximum,
                                         rdb.rule = input$rdbrule,
                                         rdbrange = c(as.numeric(input$rdbrange[1]),as.numeric(input$rdbrange[2])),
                                         ionmode = input$ionmode,
                                         all.elements = all.elements,
                                         error.ppm = input$delta.error)
      comp <- unlist(lapply(comp.object,`[[`, 1))
      theo.mass <- as.numeric(unlist(lapply(comp.object,`[[`, 2)))
      
      #extract numbers of each element
      c <- extract_num_elements("C",comp) + extract_num_elements("C[13]",comp) + extract_num_elements("C[14]",comp)
      h <- extract_num_elements("H",comp) + extract_num_elements("D[2]",comp) + extract_num_elements("T[3]",comp)
      o <- extract_num_elements("O",comp) + extract_num_elements("O[17]",comp) + extract_num_elements("O[18]",comp)
      n <- extract_num_elements("N",comp) + extract_num_elements("N[15]",comp)
      p <- extract_num_elements("P",comp)
      na <- extract_num_elements("Na",comp)
      cl <- extract_num_elements("Cl[35]",comp) + extract_num_elements("Cl[37]",comp)
    }else{
      #input is not numeric
      c <- extract_num_elements("C",input$comp)
      h <- extract_num_elements("H",input$comp)
      o <- extract_num_elements("O",input$comp)
      n <- extract_num_elements("N",input$comp)
      p <- extract_num_elements("P",input$comp)
      na <- extract_num_elements("Na",input$comp)
      cl <- extract_num_elements("Cl",input$comp)
    }
    
    str1 <- paste0("Carbon: ",c)
    str2 <- paste0("Hydrogen: ",h)
    str3 <- paste0("Oxygen: ",o)
    str4 <- paste0("Nitrogen: ",n)
    str5 <- paste0("Phosphorus: ",p)
    str6 <- paste0("Sodium: ",na)
    str7 <- paste0("Chlorine: ",cl)
    str8<-""
    str9<-""
    rdb.equiv = ( c - ((h+cl)/2) + ((n+p)/2) +1 )
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
      if(species.result %in% pls$gen.structure){
        if(input$strain != "n/a"  && species.result %in% pls[pls$strain == input$strain , colnames(pls) == "gen.structure"]){
          pls2 <- pls[pls$strain == input$strain,]
          c2 <- unname(unlist(pls2[ match(species.result,pls2$gen.structure) , 2:ncol(pls2) ]))
          c2 <- c2[c2 != ""]
          c2 <- c2[-length(c2)]
          c2 <- c(c2,paste0(" (",input$strain," specific)"))
          if(length(c2) <= 2){
            species.result <- paste0(species.result ,"</b>","; most likely exact structure: ", paste0(c2,collapse = "") )
          }else{
            species.result <- paste0(species.result ,"</b>","; most likely exact structures: ", paste(c2, collapse = ", ") )
          }
        }else{
          pls2 <- pls[pls$strain == "",]
          if(species.result %in% pls2$gen.structure){
            c2 <- unname(unlist(pls2[ match(species.result,pls2$gen.structure) , 2:ncol(pls2) ]))
            c2 <- c2[c2 != ""]
            if(length(c2) <= 1){
              species.result <- paste0(species.result ,"</b>","; most likely exact structure: ", paste0(c2,collapse = "") )
            }else{
              species.result <- paste0(species.result ,"</b>","; most likely exact structures: ", paste(c2, collapse = ", ") )
            }
          }
        }
        str12<-paste0("Species is: ","<b>",species.result)
      }else{
        if(!grepl("^L",species.result) & !grepl("FA",species.result)){
          str12<-paste0("Species is: ","<b>",species.result,"</b>",", exact structure is unknown.")
        }else{
          str12<-paste0("Species is: ","<b>",species.result)
        }
      }
    }else{
      str12<-paste0("Species not detected")
    }
    
    if(!is.na(suppressWarnings(as.numeric(input$comp)))){
      tag1 <- h1(strong("4. View results:"), style = "font-size:22.5px;")
      str1 <- paste0("Composition:")
      str2 <- ""
      str3 <- comp
      str4 <- ""
      str5 <- paste0("Theoretical mass is: ",theo.mass)
      str6 <- ""
      str7 <- paste0("Delta is: ",format(round(((theo.mass-as.numeric(input$comp))/theo.mass)*10^6,2),nsmall=2)," ppm")
    }else{
      tag1 <- h1(strong("3. View results:"), style = "font-size:22.5px;")
    }

    if(input$comp == ""){
      url2 <- a("this repository. ", href="https://github.com/briankleiboeker/SLAT", target="_blank")
      url3 <- a("the associated R package.", href = "https://github.com/briankleiboeker/slatR", target="_blank")
      tagList(
        p(HTML(paste("","This webapp was designed for use with high-resolution Orbitrap ESI-MS data with the goal of accurately assigining lipid class and structure to user-provided elemental compositions or m/z values, but should also be compatible with most any high-resolution mass spectromerty data provided in a standard format.","", sep = '<br/>'))),
        tags$img(src = "slat_logo.png", height = "300px", width = "500px",alt = "somethingwentwrong",deleteFile=FALSE),
        p(HTML(paste("","Report problems, suggestions, or other feedback to bkleiboeker [at] wustl.edu","", sep = '<br/>'))),
        p("Instructions for usage, source code, and a test dataset are available at ",url2,""," For reproducible, high-throughput assignment of lipid class and structure, check out ",url3),
        p(HTML(paste("","","Input an elemental composition or m/z value to begin",sep = '<br/>')))
      )
    }else{
      url3 <- a("this repository. ", href="https://github.com/briankleiboeker/SLAT", target="_blank")
      tagList(
      p(HTML(paste("Report problems, suggestions, or other feedback to bkleiboeker [at] wustl.edu","","Instructions for usage, source code, and a test dataset are available at ",url3,sep = '<br/>'))),
      tag1,
      p(HTML(paste(" ",str1, str2, str3, str4, str5, str6, str7,str8,str9,str10,str11,str12,sep = '<br/>')))
      )
    }
    
  })
  
  output$analyzed_data1 <- analyzed_data
  output$analyzed_data2 <- analyzed_data
  
  
  ############
  imported <- import_copypaste_server("myid",
                                      btn_show_data = F,
                                      trigger_return = "change")
  
  output$status <- renderPrint({
    if(!is.null(imported$status())){
      if(imported$status() == "success"){
        if(!("composition" %in% colnames(data.full()))){
          if(length(which(colnames(data.full()) %in% relevant.col.names)) == 0 ){
            HTML(paste0("Import successful, but no column names detected. ", "Assuming that the rightmost non-empty column contains compositions. Filter output by abs(delta) cutoff feature is not available.",ifelse( any(data.full()$structure %in% pls$gen.structure) ,""," Exact structure is unknown for all assigned structures, so exact structure columns are not being shown."),sep = '<br/>'))
          }else{
            HTML(paste0("Import successful, but composition column not detected. ", "Assuming that the rightmost non-empty column contains compositions. ",ifelse("delta" %in% colnames(data.full()),"Delta column detected, so filtering by abs(delta) cutoff feature is enabled. ","Delta column not detected, so filtering by abs(delta) cutoff feature is disabled. "),ifelse( any(data.full()$structure %in% pls$gen.structure) ,""," Exact structure is unknown for all assigned species, so exact structure columns are not being shown."),sep = '<br/>'))
          }
        }else{
          HTML(paste0("Import successful. Detected columns: ", paste( colnames(data.full())[!(colnames(data.full()) %in% c("structure","strain.specific.assignment")) & !grepl("exact.structure",colnames(data.full()))],collapse = ', '),ifelse("delta" %in% colnames(data.full()),". Delta column detected, so filtering by abs(delta) cutoff feature is enabled. ",". Delta column not detected, so filtering by abs(delta) cutoff feature is disabled. "),ifelse( any(data.full()$structure %in% pls$gen.structure) ,""," Exact structure is unknown for all assigned structures, so exact structure columns are not being shown.")))
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
      
      df$structure <- mapply(assign_species,df$c,df$h,df$o,df$n,df$p,df$na,df$cl,input$bulkionmode,rep(list(bulkadducts),nrow(df)),( df$c - ((df$h+df$cl)/2) + ((df$n+df$p)/2) +1 ),input$bulkdomain,rep(list(bulklois),nrow(df)),input$bulkmaxdblbnds,SIMPLIFY = T)
      df$structure <- ifelse(is.na(df$structure),"",df$structure)
      df<-df[, !(colnames(df) %in% c("c","h","o","n","p","na","cl"))]
      
      
      
      
      
      
      if(any(df$structure %in% pls$gen.structure)){
        
        if(input$bulkstrain != "n/a"  && any(df$structure %in% pls[pls$strain == input$bulkstrain , colnames(pls) == "gen.structure"] ) ){
          pls2 <- pls[pls$strain == input$bulkstrain,]

          df$exact.structure.1 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[1],"" )) )
          df$exact.structure.2 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[2],"" )) )
          df$exact.structure.3 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[3],"" )) )
          df$exact.structure.4 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[4],"" )) )
          df$exact.structure.5 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[5],"" )) )
          df$exact.structure.6 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[6],"" )) )
          df$exact.structure.7 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[7],"" )) )
          
          df$strain.specific.assignment <- ifelse(df$exact.structure.1 != "","Y","N")
          
          df$exact.structure.1 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pls$gen.structure,
                                                                     unname(unlist(pls[match(y,pls$gen.structure),2:ncol(pls)]))[1],
                                                                     x)},
                                                df$exact.structure.1, # this is x
                                                df$structure # this is y
                                                ))
          df$exact.structure.2 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pls$gen.structure,
                                                                     unname(unlist(pls[match(y,pls$gen.structure),2:ncol(pls)]))[2],
                                                                     x)},
                                                df$exact.structure.2, # this is x
                                                df$structure # this is y
                                                ))
          df$exact.structure.3 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pls$gen.structure,
                                                                     unname(unlist(pls[match(y,pls$gen.structure),2:ncol(pls)]))[3],
                                                                     x)},
                                                df$exact.structure.3, # this is x
                                                df$structure # this is y
                                                ))
          df$exact.structure.4 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pls$gen.structure,
                                                                     unname(unlist(pls[match(y,pls$gen.structure),2:ncol(pls)]))[4],
                                                                     x)},
                                                df$exact.structure.4, # this is x
                                                df$structure # this is y
                                                ))
          df$exact.structure.5 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pls$gen.structure,
                                                                     unname(unlist(pls[match(y,pls$gen.structure),2:ncol(pls)]))[5],
                                                                     x)},
                                                df$exact.structure.5, # this is x
                                                df$structure # this is y
                                                ))
          df$exact.structure.6 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pls$gen.structure,
                                                                     unname(unlist(pls[match(y,pls$gen.structure),2:ncol(pls)]))[6],
                                                                     x)},
                                                df$exact.structure.6, # this is x
                                                df$structure # this is y
                                                ))
          df$exact.structure.7 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pls$gen.structure,
                                                                     unname(unlist(pls[match(y,pls$gen.structure),2:ncol(pls)]))[7],
                                                                     x)},
                                                df$exact.structure.7, # this is x
                                                df$structure # this is y
                                                ))
        }else{
          pls2 <- pls[pls$strain == "",]
          df$exact.structure.1 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[1],"" )) )
          df$exact.structure.2 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[2],"" )) )
          df$exact.structure.3 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[3],"" )) )
          df$exact.structure.4 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[4],"" )) )
          df$exact.structure.5 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[5],"" )) )
          df$exact.structure.6 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[6],"" )) )
          df$exact.structure.7 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[7],"" )) )
        }
        

        if( all(df$exact.structure.7 == "") ){
          df<-df[, !(colnames(df) %in% c("exact.structure.7"))]
        }
        if( all(df$exact.structure.6 == "") ){
          df<-df[, !(colnames(df) %in% c("exact.structure.6"))]
        }
        if( all(df$exact.structure.5 == "") ){
          df<-df[, !(colnames(df) %in% c("exact.structure.5"))]
        }
        if( all(df$exact.structure.4 == "") ){
          df<-df[, !(colnames(df) %in% c("exact.structure.4"))]
        }
        if( all(df$exact.structure.3 == "") ){
          df<-df[, !(colnames(df) %in% c("exact.structure.3"))]
        }
        if( all(df$exact.structure.2 == "") ){
          df<-df[, !(colnames(df) %in% c("exact.structure.2"))]
        }
        if( all(df$exact.structure.1 == "") ){
          df<-df[, !(colnames(df) %in% c("exact.structure.1"))]
        }
        
      }
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
  
  ###############
  imported2 <- import_copypaste_server("myid2",
                                       btn_show_data = F,
                                       trigger_return = "change")
  imported3 <- import_copypaste_server("myid3",
                                       btn_show_data = F,
                                       trigger_return = "change")
  
  output$status2 <- renderPrint({
    if(!is.null(imported2$status()) | !is.null(imported3$status())){
      if(!is.null(imported2$status())){
        import.successful <- (imported2$status() == "success")
      }else{
        if(!is.null(imported3$status())){
          import.successful <- (imported3$status() == "success")
        }
      }
      if(import.successful){
        if(!("mz" %in% colnames(data.full2()))){
          if(length(which(colnames(data.full2()) %in% relevant.col.names)) == 0 ){
            HTML(paste0("Import successful, but no column names detected. ", "Assuming that the leftmost non-empty column contains m/z values.",ifelse( any(data.full2()$structure %in% pls$gen.structure) ,""," Exact structure is unknown for all assigned structures, so exact structure columns are not being shown."),sep = '<br/>'))
          }else{
            HTML(paste0("Import successful, but m/z column not detected. ", "Assuming that the leftmost non-empty column contains m/z values. ",ifelse( any(data.full2()$structure %in% pls$gen.structure) ,""," Exact structure is unknown for all assigned species, so exact structure columns are not being shown."),sep = '<br/>'))
          }
        }else{
          HTML(paste0("Import successful. Detected columns: ", paste( colnames(data.full2())[ !(colnames(data.full2()) %in% c("structure","composition","theo.mass","rdb","delta","strain.specific.assignment")) & !grepl("exact.structure",colnames(data.full2()))],collapse = ', '),ifelse( any(data.full2()$structure %in% pls$gen.structure) ,"",". Exact structure is unknown for all assigned structures, so exact structure columns are not being shown.")))
        }
      }else{
        if(!is.null(imported2$status())){
          imported2$status()
        }else{
          imported3$status()
        }
        
      }
    }else{
      imported2$status()
    }
  })
  
  output$status3 <- renderPrint({
    if(!is.null(imported3$status())){
      if(imported3$status() == "success"){
        imported3$status()
        if(!("mz" %in% colnames(data.full3()))){
          if(length(which(colnames(data.full3()) %in% relevant.col.names)) == 0 ){
            HTML(paste0("Import successful, but no column names detected. ", "Assuming that the leftmost non-empty column contains m/z values.",sep = '<br/>'))
          }else{
            HTML(paste0("Import successful, but m/z column not detected. ", "Assuming that the leftmost non-empty column contains m/z values. ",sep = '<br/>'))
          }
        }else{
          HTML(paste0("Import successful. Detected columns: ", paste( colnames(data.full3())[ !(colnames(data.full3()) %in% c("corr.mz","intensity.plot.col"))],collapse = ', ')))
        }
      }else{
        imported3$status()
      }
    }else{
      imported3$status()
    }
  })
  
  output$name2 <- renderPrint({
    imported2$name()
  })
  
  output$name3 <- renderPrint({
    imported3$name()
  })
  
  ############
  
  alignment_table <- reactiveVal(data.frame(mz = numeric(),
                                            theo.mass = numeric()))
  
  output$alignment_shiny_table <- renderDT({
      datatable(alignment_table(), selection = 'single', options = list(dom = 't'),editable = T)
  })

  observeEvent(input$alignment_delete_btn, {
    t = alignment_table()
    if (!is.null(input$alignment_shiny_table_rows_selected)) {
      t <- t[-as.numeric(input$alignment_shiny_table_rows_selected),]
    }
    alignment_table(t)
  })
  
  observeEvent(input$alignment_shiny_table_cell_edit, {
    at <- alignment_table()
    at[input$alignment_shiny_table_cell_edit$row,input$alignment_shiny_table_cell_edit$col] <- input$alignment_shiny_table_cell_edit$value
    alignment_table(at)
  })
  

  observeEvent(event_data("plotly_click", source = "plot1"), {
    d <- event_data("plotly_click", source = "plot1")
    add_row <- data.frame(mz = d$x,
                          theo.mass = NA)
    if( !(add_row$mz %in% alignment_table()$mz) ){
      alignment_table(rbind(add_row, alignment_table()))
    }
  })
  
  
  output$plot1 <- renderPlotly({
    if(!is.null(tab3data())){
      df <- data.full3()
      df$mz <- round(df$mz, digits = 8)
      divider <- 100
       p1 <- ggplot(df, aes(x = mz, y = intensity.plot.col))+
          geom_point()+
          theme_bw()+
          ylab("Intensity")+
          xlab("m/z")+
          annotate("text",
                   #data = dplyr::slice_head(dplyr::arrange(df,desc(intensity.plot.col)), n = max(10,floor(nrow(df)/100))),
                   x = dplyr::slice_head(dplyr::arrange(df,desc(intensity.plot.col)), n = max(10,floor(nrow(df)/divider)))$mz,
                   y = dplyr::slice_head(dplyr::arrange(df,desc(intensity.plot.col)), n = max(10,floor(nrow(df)/divider)))$intensity.plot.col + max(abs(df$intensity.plot.col))/30,
                   label = dplyr::slice_head(dplyr::arrange(df,desc(intensity.plot.col)), n = max(10,floor(nrow(df)/divider)))$mz,
                   # vjust = max(abs(df$intensity.plot.col))/10,
                   size = 3
                   )+
          # geom_text(data = dplyr::slice_head(dplyr::arrange(df,desc(intensity.plot.col)), n = max(10,floor(nrow(df)/100)) ),
          #                          aes(label = mz,x = mz, y = intensity.plot.col + max(abs(df$intensity.plot.col))/50))+
          theme(text=element_text(size=15))
       p1 %>% ggplotly(source = "plot1") %>% 
         event_register('plotly_click') %>% 
         config(displayModeBar = FALSE)

    }
  })
  
  
  data.full3 <- reactive({
    if(!is.null(tab3data())){
      df <- tab3data()
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

      if("mz" %in% colnames(df)){
        df <- df[df$mz != "" & !is.na(df$mz) & !is.null(df$mz),]
      }else{
        df <- df[df[,1] != "" & !is.na(df[,1]) & !is.null(df[,1]),]
        colnames(df)[1] <- "mz"
      }
      if( !("rel.intensity" %in% colnames(df)) && !("intensity" %in% colnames(df))){
        colnames(df)[2]<-"intensity"
        df$intensity.plot.col <- df$intensity
      }else{
        if("rel.intensity" %in% colnames(df)){
          df$intensity.plot.col <- df$rel.intensity
        }else{
          df$intensity.plot.col <- df$intensity
        }
      }
        
      if(nrow(alignment_table()) > 0 && all(!is.na(alignment_table()$theo.mass)) ){
        if(nrow(alignment_table()) == 1){
          df$corr.mz <- as.numeric(df$mz) + (as.numeric(alignment_table()$theo.mass[1]) - as.numeric(alignment_table()$mz[1]))
        }else{
          yvals <- (as.numeric(alignment_table()$theo.mass) - as.numeric(alignment_table()$mz))
          xvals <- as.numeric(alignment_table()$mz)
          s <- summary(lm(yvals ~ xvals))$coefficients
          yint <- as.numeric(s[rownames(s) == "(Intercept)",colnames(s) == "Estimate"])
          slope <- as.numeric(s[rownames(s) == "xvals",colnames(s) == "Estimate"])
          df$corr.mz <- as.numeric(df$mz) + ((slope*as.numeric(df$mz))+yint)
        }
      }else{
        df$corr.mz <- rep(NA,nrow(df))
      }
      df
    }
    else{
      tab3data()
    }
    
  })
  
  
  tab2data <- reactive({
      return(imported2$data())
  })
  
  tab3data <- reactive({
      return(imported3$data())
  })
  
  data.full2 <- reactive({
    transferdata_value <- ifelse(is.null(input$transferdata),
                                 F,
                                 input$transferdata == "y")
    if( !is.null( tab2data() ) | !is.null(tab3data()) ){
      if(transferdata_value | is.null(tab2data() )){
        df <- tab3data()
        if(!is.null(df)){
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
          
          if("mz" %in% colnames(df)){
            df <- df[df$mz != "" & !is.na(df$mz) & !is.null(df$mz),]
          }else{
            df <- df[df[,1] != "" & !is.na(df[,1]) & !is.null(df[,1]),]
            colnames(df)[1] <- "mz"
          }
          
          if(nrow(alignment_table()) > 0 && all(!is.na(alignment_table()$theo.mass)) ){
            if(nrow(alignment_table()) == 1){
              df$mz <- as.numeric(df$mz) + (as.numeric(alignment_table()$theo.mass[1]) - as.numeric(alignment_table()$mz[1]))
            }else{
              yvals <- (as.numeric(alignment_table()$theo.mass) - as.numeric(alignment_table()$mz))
              xvals <- as.numeric(alignment_table()$mz)
              s <- summary(lm(yvals ~ xvals))$coefficients
              yint <- as.numeric(s[rownames(s) == "(Intercept)",colnames(s) == "Estimate"])
              slope <- as.numeric(s[rownames(s) == "xvals",colnames(s) == "Estimate"])
              df$mz <- as.numeric(df$mz) + ((slope*as.numeric(df$mz))+yint)
            }
          }
        }
        
      }else{
        df <- tab2data()
      }

      
      withProgress(message = "Assigning composition to m/z values",value = 0, {
      while( most(is.na(df[,ncol(df)])) | most(df[,ncol(df)] == "") | most(is.null(df[,ncol(df)]))){
        df<-df[,1:(ncol(df)-1)]
      }
      while(most(is.na(df[,1])) | most(df[,1] == "") | most(is.null(df[,1]))){
        df<-df[,2:ncol(df)]
      }
      
      colnames(df)[(unlist(lapply(colnames(df),grepl,pattern = "m",ignore.case = T,perl = T)) & unlist(lapply(colnames(df),grepl,pattern = "z",ignore.case = T,perl = T)))]<-"mz"
      colnames(df)[(unlist(lapply(colnames(df),grepl,pattern = "int",ignore.case = T,perl = T)) & !unlist(lapply(colnames(df),grepl,pattern = "rel",ignore.case = T,perl = T)))]<-"intensity"
      colnames(df)[unlist(lapply(colnames(df),grepl,pattern = "rel",ignore.case = T,perl = T))]<-"rel.intensity"
      colnames(df)[unlist(lapply(colnames(df),grepl,pattern = "theo",ignore.case = T,perl = T))]<-"theo.mass"
      colnames(df)[( unlist(lapply(colnames(df),grepl,pattern = "delta",ignore.case = T,perl = T)) | unlist(lapply(colnames(df),grepl,pattern = "mmu",ignore.case = T,perl = T)) | unlist(lapply(colnames(df),grepl,pattern = "amu",ignore.case = T,perl = T)) | unlist(lapply(colnames(df),grepl,pattern = "ppm",ignore.case = T,perl = T)) )]<-"delta"
      colnames(df)[unlist(lapply(colnames(df),grepl,pattern = "rdb",ignore.case = T,perl = T))]<-"rdb.equiv"
      colnames(df)[unlist(lapply(colnames(df),grepl,pattern = "comp",ignore.case = T,perl = T))]<-"composition"
      
      if("mz" %in% colnames(df)){
        df <- df[df$mz != "" & !is.na(df$mz) & !is.null(df$mz),]
        comp.object <- assign_comp_from_mz(mz = df$mz,
                                           elements = this_table2()$element,
                                           isotope = this_table2()$isotope,
                                           atomic.mass = this_table2()$atomic.mass,
                                           mins.list = this_table2()$minimum,
                                           maxs.list = this_table2()$maximum,
                                           rdb.rule = input$bulkrdbrule,
                                           rdbrange = c(as.numeric(input$bulkrdbrange[1]),as.numeric(input$bulkrdbrange[2])),
                                           ionmode = input$bulkionmode,
                                           all.elements = all.elements,
                                           error.ppm = input$bulk.delta.error)
      }else{
        df <- df[df[,1] != "" & !is.na(df[,1]) & !is.null(df[,1]),]
        comp.object <- assign_comp_from_mz(mz = df[,1],
                                           elements = this_table2()$element,
                                           isotope = this_table2()$isotope,
                                           atomic.mass = this_table2()$atomic.mass,
                                           mins.list = this_table2()$minimum,
                                           maxs.list = this_table2()$maximum,
                                           rdb.rule = input$bulkrdbrule,
                                           rdbrange = c(as.numeric(input$bulkrdbrange[1]),as.numeric(input$bulkrdbrange[2])),
                                           ionmode = input$bulkionmode,
                                           all.elements = all.elements,
                                           error.ppm = input$bulk.delta.error)
      }
      
      df$composition <- unlist(lapply(comp.object,`[[`, 1))
      df$theo.mass <- as.numeric(unlist(lapply(comp.object,`[[`, 2)))
      df$theo.mass[is.na(df$theo.mass)] <- ""
      df$delta <- as.numeric(unlist(lapply(comp.object,`[[`, 3)))
      df$delta[is.na(df$delta)] <- ""
      df$c <- return_num_elements("C",df$composition) + return_num_elements("C[13]",df$composition) + return_num_elements("C[14]",df$composition)
      df$h <- return_num_elements("H",df$composition) + return_num_elements("D[2]",df$composition) + return_num_elements("T[3]",df$composition)
      df$o <- return_num_elements("O",df$composition) + return_num_elements("O[17]",df$composition) + return_num_elements("O[18]",df$composition)
      df$n <- return_num_elements("N",df$composition) + return_num_elements("N[15]",df$composition)
      df$p <- return_num_elements("P",df$composition)
      df$na <- return_num_elements("Na",df$composition)
      df$cl <- return_num_elements("Cl[35]",df$composition) + return_num_elements("Cl[37]",df$composition)
      df$rdb <- ifelse(df$theo.mass != "",df$c - ((df$h + df$cl) / 2) + ((df$n+df$p)/2) + 1,"")

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

      df$structure <- mapply(assign_species,df$c,df$h,df$o,df$n,df$p,df$na,df$cl,input$bulkionmode,rep(list(bulkadducts),nrow(df)),( df$c - ((df$h + df$cl)/2) + ((df$n+df$p)/2) +1 ),input$bulkdomain,rep(list(bulklois),nrow(df)),input$bulkmaxdblbnds,SIMPLIFY = T)
      df$structure <- ifelse(is.na(df$structure),"",df$structure)
      df<-df[, !(colnames(df) %in% c("c","h","o","n","p","na","cl"))]

      if(any(df$structure %in% pls$gen.structure)){

        if(input$bulkstrain != "n/a"  && any(df$structure %in% pls[pls$strain == input$bulkstrain , colnames(pls) == "gen.structure"] ) ){
          pls2 <- pls[pls$strain == input$bulkstrain,]
          
          df$exact.structure.1 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[1],"" )) )
          df$exact.structure.2 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[2],"" )) )
          df$exact.structure.3 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[3],"" )) )
          df$exact.structure.4 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[4],"" )) )
          df$exact.structure.5 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[5],"" )) )
          df$exact.structure.6 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[6],"" )) )
          df$exact.structure.7 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[7],"" )) )
          
          df$strain.specific.assignment <- ifelse(df$exact.structure.1 != "","Y","N")
          
          df$exact.structure.1 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pls$gen.structure,
                                                                     unname(unlist(pls[match(y,pls$gen.structure),2:ncol(pls)]))[1],
                                                                     x)},
                                                df$exact.structure.1, # this is x
                                                df$structure # this is y
          ))
          df$exact.structure.2 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pls$gen.structure,
                                                                     unname(unlist(pls[match(y,pls$gen.structure),2:ncol(pls)]))[2],
                                                                     x)},
                                                df$exact.structure.2, # this is x
                                                df$structure # this is y
          ))
          df$exact.structure.3 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pls$gen.structure,
                                                                     unname(unlist(pls[match(y,pls$gen.structure),2:ncol(pls)]))[3],
                                                                     x)},
                                                df$exact.structure.3, # this is x
                                                df$structure # this is y
          ))
          df$exact.structure.4 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pls$gen.structure,
                                                                     unname(unlist(pls[match(y,pls$gen.structure),2:ncol(pls)]))[4],
                                                                     x)},
                                                df$exact.structure.4, # this is x
                                                df$structure # this is y
          ))
          df$exact.structure.5 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pls$gen.structure,
                                                                     unname(unlist(pls[match(y,pls$gen.structure),2:ncol(pls)]))[5],
                                                                     x)},
                                                df$exact.structure.5, # this is x
                                                df$structure # this is y
          ))
          df$exact.structure.6 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pls$gen.structure,
                                                                     unname(unlist(pls[match(y,pls$gen.structure),2:ncol(pls)]))[6],
                                                                     x)},
                                                df$exact.structure.6, # this is x
                                                df$structure # this is y
          ))
          df$exact.structure.7 <- unlist(mapply(function(x,y){ifelse(x == "" && y %in% pls$gen.structure,
                                                                     unname(unlist(pls[match(y,pls$gen.structure),2:ncol(pls)]))[7],
                                                                     x)},
                                                df$exact.structure.7, # this is x
                                                df$structure # this is y
          ))
        }else{
          pls2 <- pls[pls$strain == "",]
          df$exact.structure.1 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[1],"" )) )
          df$exact.structure.2 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[2],"" )) )
          df$exact.structure.3 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[3],"" )) )
          df$exact.structure.4 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[4],"" )) )
          df$exact.structure.5 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[5],"" )) )
          df$exact.structure.6 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[6],"" )) )
          df$exact.structure.7 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[7],"" )) )
        }
        

        if( all(df$exact.structure.7 == "") ){
          df<-df[, !(colnames(df) %in% c("exact.structure.7"))]
        }
        if( all(df$exact.structure.6 == "") ){
          df<-df[, !(colnames(df) %in% c("exact.structure.6"))]
        }
        if( all(df$exact.structure.5 == "") ){
          df<-df[, !(colnames(df) %in% c("exact.structure.5"))]
        }
        if( all(df$exact.structure.4 == "") ){
          df<-df[, !(colnames(df) %in% c("exact.structure.4"))]
        }
        if( all(df$exact.structure.3 == "") ){
          df<-df[, !(colnames(df) %in% c("exact.structure.3"))]
        }
        if( all(df$exact.structure.2 == "") ){
          df<-df[, !(colnames(df) %in% c("exact.structure.2"))]
        }
        if( all(df$exact.structure.1 == "") ){
          df<-df[, !(colnames(df) %in% c("exact.structure.1"))]
        }
        }
       })
       
      df
      
    }
    else{
      tab2data()
    }
  }) %>% bindEvent(input$bulk.analyze.button)
  
  
  output$data2 <- renderPrint({
    if(is.null(tab2data()) & is.null(tab3data())){
      tab2data()
    }else{
      if(input$returnAll == "returnassignedpeaks"){
        data.full2()[data.full2()$structure != "",] |> head(n = 10)
        }else{
          data.full2() |> head(n = 10)
          }
    }
      
  })
  
  output$clip <- renderUI({
    if(!is.null(data.full3())){
      rclipButton(inputId = "clipbtn",
                  label = "Copy aligned data to clipboard",
                  clipText = readr::format_tsv(suppressWarnings(dplyr::rename(dplyr::select(data.full3()[,!(colnames(data.full3()) %in% c("intensity.plot.col"))],!c("mz")),"mz" = "corr.mz") )),
                  icon = icon("clipboard"))
      
    }else{
      rclipButton(inputId = "clipbtn",
                  label = "Copy aligned data to clipboard", 
                  clipText = "something went wrong",
                  icon = icon("clipboard"))
      
    }

  })
  
  output$download2 <- downloadHandler(
    filename = function(){
      paste0(imported2$name(),".csv")
    },
    content = function(file){
      if(( !is.null(tab2data()) | !is.null(tab3data()) ) && input$returnAll == "returnassignedpeaks"){
        write.csv(data.full2()[data.full2()$structure != "",],file,row.names = F)
      }else{
        write.csv(data.full2(),file,row.names = F)
      }
    }
  )
  
  output$download3 <- downloadHandler(
    filename = function(){
      paste0(imported3$name(),".csv")
    },
    content = function(file){
      write.csv(data.full3()[,!(colnames(data.full3()) %in% c("intensity.plot.col"))],file,row.names = F)
    }
  )
  
}

shinyApp(ui = ui, server = server)
