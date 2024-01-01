library(shiny)
library(datamods)
library(DT)
library(ggplot2)
library(rclipboard)
library(plotly)

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

assign_species <- function(c,h,o,n,p,s,na,cl,li,ion.mode,adducts,rdb,domain,lois,max.dbl.bnds,strain){
  if(is.na(rdb) | is.null(adducts)){
    return(NA)
  }
  if(ion.mode == "neutral"){
    return(structure_from_exact_comp(c,h,o,n,p,s,c("nonacylglycerols","acylglycerols"),domain,lois,as.numeric(max.dbl.bnds),strain))
  }
  if(ion.mode == "pos.ion"){
    if("m.plus.sodium" %in% adducts){
      if(na > 0){
        return(structure_from_exact_comp(c,h,o,n,p,s,c("nonacylglycerols","acylglycerols","galactosyldiacylglycerols"),domain,lois,as.numeric(max.dbl.bnds),strain))
        # look for sodium adduct species -- ALL
        # assign exact composition by changing nothing
        # regardless if find one or not, return the result here
      }
      # else, keep going
    }
    if("m.plus.lithium" %in% adducts){
      if(li > 0){
        return(structure_from_exact_comp(c,h,o,n,p,s,c("nonacylglycerols","acylglycerols","galactosyldiacylglycerols"),domain,lois,as.numeric(max.dbl.bnds),strain))
        # look for lithium adduct species -- ALL
        # assign exact composition by changing nothing
        # regardless if find one or not, return the result here
      }
      # else, keep going
    }
    if("m.plus.ammonia" %in% adducts){
      if(n > 0){
        result <- structure_from_exact_comp(c,h-4,o,n-1,p,s,c("acylglycerols"),domain,lois,as.numeric(max.dbl.bnds),strain)
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
      return(structure_from_exact_comp(c,h-1,o,n,p,s,c("nonacylglycerols"),domain,lois,as.numeric(max.dbl.bnds),strain))
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
        return(structure_from_exact_comp(c,h,o,n,p,s,c("nonacylglycerols","galactosyldiacylglycerols"),domain,lois,as.numeric(max.dbl.bnds),strain))
        # look for chloride adduct species
        # assign exact composition by changing nothing
        # regardless if find one or not, return the result here
      }
      # else, keep going
    }
    if("m.minus.h" %in% adducts){
      result <- structure_from_exact_comp(c,h+1,o,n,p,s,c("nonacylglycerols","galactosyldiacylglycerols"),domain,lois,as.numeric(max.dbl.bnds),strain)
      if(!is.na(result)){
        return(result)
      }
      # look for lipid matches
      # do this by adding 1H and seeing if the species could exist
      # if it does exist return m.minus.h
      # if it doesn't exist -> return some clever paste0 showing what we searched for, beginning with 'Species not detected (traceback:'
    }
    if("m.minus.2h" %in% adducts){
      result <- structure_from_exact_comp(c,h+2,o,n,p,s,c("nonacylglycerols"),domain,lois,as.numeric(max.dbl.bnds),strain)
      if(!is.na(result)){
        return(result)
      }
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
structure_from_exact_comp <- function(c,h,o,n,p,s,which_to_look_for,domain,lois,max.dbl.bnds,strain){
  if(domain == "bact"){
    lyso.carbon.cutoff.value <- 19
  }else{
    lyso.carbon.cutoff.value <- 27
  }
  rdb.equiv = ( c - (h/2) + ((n+p)/2) +1 )
  if(rdb.equiv %% 1 != 0 ){
    return(NA)
  }
  
  if(domain == "euk" | strain == "e.coli" | strain == "listeria" | strain == "streptococcus" | strain == "n/a"){
    if( "acylglycerols" %in% which_to_look_for){
      if(n == 0 && p == 0 && o == 5 && c >= 19 && s == 0){
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
      if(n == 0 && p == 0 && o == 6 && c >= 27 && s == 0){
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
    if( ("galactosyldiacylglycerols" %in% which_to_look_for)  &  (strain == "listeria" | strain == "streptococcus") ){
      if(n == 0 && p == 0 && c > 29 && (o == 10 | o == 15)  && s == 0){
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
      if(n == 0 && s == 0){
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
      if(n >= 1 && s == 0 && ("galcer" %in% lois) ){
        #look for GalCer and LacCer, but do not ever return(NA) if they are not found
        if( (o == 8 | o == 9) && rdb.equiv >= 2){
          #GalCer
          if(o == 8){
            n.double.bonds <-  (rdb.equiv - 2)
            n.fa.carbons <- (c - 6)
            if(n.double.bonds < max.dbl.bnds){
              return(paste0("GalCer(",n.fa.carbons,":",n.double.bonds,")"))
            }
          }else{
            n.double.bonds <-  (rdb.equiv - 2)
            n.fa.carbons <- (c - 6)
            if(n.double.bonds < max.dbl.bnds){
              return(paste0("GalCer(h",n.fa.carbons,":",n.double.bonds,")"))
            }
          }

        }
        if(o == 13 && rdb.equiv >= 3){
          #LacCer
          n.double.bonds <-  (rdb.equiv - 3)
          n.fa.carbons <- (c - 12)
          if(n.double.bonds < max.dbl.bnds){
            return(paste0("LacCer(",n.fa.carbons,":",n.double.bonds,")"))
          }
        }
      }
      if(n >= 1 && s == 0 && ("gng" %in% lois) ){
        #look for Gangliosides, but do not ever return(NA) if they are not found
        if(o == 18 && rdb.equiv >= 4){
          #GB3; GA2
          if(rdb.equiv <= 5){
            #most likely a GB3
            n.double.bonds <-  (rdb.equiv - 4)
            n.fa.carbons <- (c - 18)
            if(n.double.bonds < max.dbl.bnds){
              return(paste0("GB3(",n.fa.carbons,":",n.double.bonds,")"))
            }
          }else{
            #could be either, but best to assume it is a GA2
            n.double.bonds <-  (rdb.equiv - 5)
            n.fa.carbons <- (c - 20)
            if(n.double.bonds < max.dbl.bnds){
              return(paste0("GA2(",n.fa.carbons,":",n.double.bonds,")"))
            }
          }
        }
        if((o == 23 | o == 24) && rdb.equiv >= 6){
          #GB4; GA1
          if(o == 23){
            n.double.bonds <-  (rdb.equiv - 6)
            n.fa.carbons <- (c - 26)
            if(n.double.bonds < max.dbl.bnds){
              return(paste0("GB4(",n.fa.carbons,":",n.double.bonds,")"))
            }
          }else{
            n.double.bonds <-  (rdb.equiv - 6)
            n.fa.carbons <- (c - 26)
            if(n.double.bonds < max.dbl.bnds){
              return(paste0("GB4(h",n.fa.carbons,":",n.double.bonds,")"))
            }
          }

        }
        if(o == 21 && rdb.equiv >= 6){
          #GM3
          n.double.bonds <-  (rdb.equiv - 6)
          n.fa.carbons <- (c - 23)
          if(n.double.bonds < max.dbl.bnds){
            return(paste0("GM3(",n.fa.carbons,":",n.double.bonds,")"))
          }
        }
        if(o == 26 && rdb.equiv >= 8){
          #GM2
          n.double.bonds <-  (rdb.equiv - 8)
          n.fa.carbons <- (c - 31)
          if(n.double.bonds < max.dbl.bnds){
            return(paste0("GM2(",n.fa.carbons,":",n.double.bonds,")"))
          }
        }
        if(o == 31 && rdb.equiv >= 9){
          #GM1
          n.double.bonds <-  (rdb.equiv - 9)
          n.fa.carbons <- (c - 37)
          if(n.double.bonds < max.dbl.bnds){
            return(paste0("GM1(",n.fa.carbons,":",n.double.bonds,")"))
          }
        }
        if(o == 27 && rdb.equiv >= 7){
          #GD3
          n.double.bonds <-  (rdb.equiv - 7)
          n.fa.carbons <- (c - 30)
          if(n.double.bonds < max.dbl.bnds){
            return(paste0("GD3(",n.fa.carbons,":",n.double.bonds,")"))
          }
        }
        if(o == 32 && rdb.equiv >= 9){
          #GD2
          n.double.bonds <-  (rdb.equiv - 9)
          n.fa.carbons <- (c - 38)
          if(n.double.bonds < max.dbl.bnds){
            return(paste0("GD2(",n.fa.carbons,":",n.double.bonds,")"))
          }
        }
        if(o == 37 && rdb.equiv >= 10){
          #GD1
          n.double.bonds <-  (rdb.equiv - 10)
          n.fa.carbons <- (c - 44)
          if(n.double.bonds < max.dbl.bnds){
            return(paste0("GD1(",n.fa.carbons,":",n.double.bonds,")"))
          }
        }
        if(o == 33 && rdb.equiv >= 8){
          #GT3
          n.double.bonds <-  (rdb.equiv - 8)
          n.fa.carbons <- (c - 37)
          if(n.double.bonds < max.dbl.bnds){
            return(paste0("GT3(",n.fa.carbons,":",n.double.bonds,")"))
          }
        }
        if(o == 38 && rdb.equiv >= 10){
          #GT2
          n.double.bonds <-  (rdb.equiv - 10)
          n.fa.carbons <- (c - 45)
          if(n.double.bonds < max.dbl.bnds){
            return(paste0("GT2(",n.fa.carbons,":",n.double.bonds,")"))
          }
        }
        if(o == 43 && rdb.equiv >= 11){
          #GT1
          n.double.bonds <-  (rdb.equiv - 11)
          n.fa.carbons <- (c - 51)
          if(n.double.bonds < max.dbl.bnds){
            return(paste0("GT1(",n.fa.carbons,":",n.double.bonds,")"))
          }
        }
        if(o == 48 && rdb.equiv >= 13){
          #GQ1
          n.double.bonds <-  (rdb.equiv - 13)
          n.fa.carbons <- (c - 59)
          if(n.double.bonds < max.dbl.bnds){
            return(paste0("GQ1(",n.fa.carbons,":",n.double.bonds,")"))
          }
        }
      }
      if(n == 2 && s == 0){
        # if(p == 1 && o == 7 && ("epc" %in% lois) ){
        #   #an epc!
        #   #make sure to keep this BEFORE the SM statements, and do not put unnecessary return(NA) statements here
        #   if(rdb.equiv >= 1){
        #     n.double.bonds <-  (rdb.equiv - 1)
        #     n.fa.carbons <- (c - 2)
        #     if(n.double.bonds < max.dbl.bnds){
        #       return(paste0("EPC(h",n.fa.carbons,":",n.double.bonds,")"))
        #     }
        #   }
        # }
        if(p == 1 && o == 6){
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
        # if(p == 1 && o == 9){
        #   #a spc!
        #   if( (rdb.equiv < 2) | !("spc" %in% lois) ){
        #     return(NA)
        #   }
        #   n.double.bonds <-  (rdb.equiv - 2)
        #   n.fa.carbons <- (c - 3)
        #   if(n.double.bonds > max.dbl.bnds){
        #     return(NA)
        #   }
        #   return(paste0("SPC(h",n.fa.carbons,":",n.double.bonds,")"))
        # }
        return(NA)
      }
      if(n == 1){
        if(p == 0 && o >= 3 && o <= 5 && s == 0){
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
        if(p == 0 && (o == 11 | o == 12) && s == 1){
          ##### sulfatides #####
          
          if(rdb.equiv < 2 | !("sulf" %in% lois)){
            return(NA)
          }
          n.double.bonds <-  (rdb.equiv - 2)
          n.fa.carbons <- (c - 6)
          if(n.double.bonds > max.dbl.bnds){
            return(NA)
          }
          if(o == 11){
            return(paste0("Sulf(",n.fa.carbons,":",n.double.bonds,")"))
          }else{
            # o == 12 here
            return(paste0("Sulf(h",n.fa.carbons,":",n.double.bonds,")"))
          }

        }
        if(p == 1 && s == 0){
          # if(o == 9 && ("gpc" %in% lois)){
          #   #an GPC!
          #   #make sure to keep this BEFORE the PS statements, and do not put return(NA) statements here
          #   #have no idea if any of this is right for RDB cutoff or n.fa.carbon cutoff
          #   if(rdb.equiv >= 1){
          #     n.double.bonds <-  (rdb.equiv - 1)
          #     n.fa.carbons <- (c - 3)
          #     if(n.double.bonds < max.dbl.bnds){
          #       return(paste0("GPC(h",n.fa.carbons,":",n.double.bonds,")"))
          #     }
          #   }
          #   
          # }
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
          # if(o == 7 && ("cer1p" %in% lois)){
          #   #a Cer-1-P
          #   #make sure to keep this BEFORE the PC/PE statements, and do not put unnecessary return(NA) statements here
          #   if(rdb.equiv >= 1){
          #     n.double.bonds <-  (rdb.equiv - 1)
          #     n.fa.carbons <- (c - 0)
          #     if(n.double.bonds < max.dbl.bnds){
          #       return(paste0("Cer-1-P(h",n.fa.carbons,":",n.double.bonds,")"))
          #     }
          #   }
          #   
          # }
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
          # if(o == 12){
          #   #an IPC!
          #   if( (rdb.equiv < 2) | !("ipc" %in% lois) ){
          #     return(NA)
          #   }
          #   n.double.bonds <-  (rdb.equiv - 2)
          #   n.fa.carbons <- (c - 6)
          #   if(n.double.bonds > max.dbl.bnds){
          #     return(NA)
          #   }
          #   return(paste0("IPC(h",n.fa.carbons,":",n.double.bonds,")"))
          # }
        }
        if(length(which_to_look_for)==1){
          return(NA)
        }
      }
    }
  }
  if(domain == "bact"){
    if(strain == "bacteroidetes" | strain == "p.gingivalis"){
      if("nonacylglycerols" %in% which_to_look_for){ #look for everything but TAG and DAG
        if(n == 0 && s == 0){
          if(p == 1){
            #pa, pg, or pi
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
          return(NA)
        }
        if(n == 2 && s == 0){
          if(p == 1 && (o == 7 | o == 8) && ("epc" %in% lois) ){
            #an epc!
            #make sure to keep this BEFORE the SM statements, and do not put unnecessary return(NA) statements here
            if(rdb.equiv >= 1 && o == 7){
              n.double.bonds <-  (rdb.equiv - 1)
              n.fa.carbons <- (c - 2)
              if(n.double.bonds < max.dbl.bnds){
                return(paste0("EPC(h",n.fa.carbons,":",n.double.bonds,")"))
              }
            }
            if(rdb.equiv >= 2 && o == 8 && strain == "p.gingivalis"){
              n.double.bonds <-  (rdb.equiv - 2)
              n.fa.carbons <- (c - 1)
              if(n.double.bonds < max.dbl.bnds){
                return(paste0("acyl-EPC(h",n.fa.carbons,":",n.double.bonds,")"))
              }
            }
          }
          if(p == 2){
            if(o >= 15 && o <= 16 && ("cer" %in% lois) ){
              #cer-PGP-cer
              if( (rdb.equiv < 2 && o == 15) | (rdb.equiv < 3 && o == 16) ){
                return(NA)
              }
              
              if(o == 15){
                n.double.bonds <-  (rdb.equiv - 2)
                n.fa.carbons <- (c - 3)
                if(n.double.bonds > max.dbl.bnds){
                  return(NA)
                }
                return(paste0("Cer-PGP-Cer(h2,",n.fa.carbons,":",n.double.bonds,")"))
              }else{
                n.double.bonds <-  (rdb.equiv - 3)
                n.fa.carbons <- (c - 3)
                if(n.double.bonds > max.dbl.bnds){
                  return(NA)
                }
                return(paste0("acyl-Cer-PGP-Cer(h2,",n.fa.carbons,":",n.double.bonds,")"))
              }
            }
            if(o == 18){
              #a dhc-pip-dhc
              if( (rdb.equiv < 3) | !("dpd" %in% lois) ){
                return(NA)
              }
              n.double.bonds <-  (rdb.equiv - 3)
              n.fa.carbons <- (c - 6)
              if(n.double.bonds > max.dbl.bnds){
                return(NA)
              }
              return(paste0("DHC-PIP-DHC(h2,",n.fa.carbons,":",n.double.bonds,")"))
            }
          }
          if(p == 0 && (o == 7 | o == 6)){
            #a GS-lipid!
            if((o == 6 && rdb.equiv < 3)  | (o == 7 && rdb.equiv < 4) | !("gslipid" %in% lois) ){
              return(NA)
            }
            
            if(o == 6){
              n.double.bonds <-  (rdb.equiv - 3)
              n.fa.carbons <- (c - 5)
              if(n.double.bonds > max.dbl.bnds){
                return(NA)
              }
              return(paste0("GS(",n.fa.carbons,":",n.double.bonds,")"))
              
            }else{
              n.double.bonds <-  (rdb.equiv - 4)
              n.fa.carbons <- (c - 5)
              if(n.double.bonds > max.dbl.bnds){
                return(NA)
              }
              return(paste0("acyl-GS(",n.fa.carbons,":",n.double.bonds,")"))
              
            }

          }
          if(p == 1 && (o == 14 | o == 13)){
            #a GS-PA!
            if((rdb.equiv < 5 && o == 13) | (rdb.equiv < 6 && o == 14) | !("gspa" %in% lois) ){
              return(NA)
            }
            if(o == 13){
              n.double.bonds <-  (rdb.equiv - 5)
              n.fa.carbons <- (c - 8)
              if(n.double.bonds > max.dbl.bnds){
                return(NA)
              }
              return(paste0("GS-PA(h",n.fa.carbons,":",n.double.bonds,")"))
              
            }else{ #o == 14, an acyl-GS-PA
              n.double.bonds <-  (rdb.equiv - 6)
              n.fa.carbons <- (c - 8)
              if(n.double.bonds > max.dbl.bnds){
                return(NA)
              }
              return(paste0("acyl-GS-PA(h",n.fa.carbons,":",n.double.bonds,")"))
            }
          }
          if(p == 1 && o == 9){
            #a spc!
            if( (rdb.equiv < 2) | !("spc" %in% lois) ){
              return(NA)
            }
            n.double.bonds <-  (rdb.equiv - 2)
            n.fa.carbons <- (c - 3)
            if(n.double.bonds > max.dbl.bnds){
              return(NA)
            }
            return(paste0("SPC(h",n.fa.carbons,":",n.double.bonds,")"))
          }
          if(p == 1 && o == 12){
            #acyl-GS-PG!
            if( (rdb.equiv < 4) | !("gspg" %in% lois) ){
              return(NA)
            }
            n.double.bonds <-  (rdb.equiv - 4)
            n.fa.carbons <- (c - 8)
            if(n.double.bonds > max.dbl.bnds){
              return(NA)
            }
            return(paste0("acyl-GS-PG(h",n.fa.carbons,":",n.double.bonds,")"))
          }
          return(NA)
        }
        if(n == 1 && s == 0){
          if(p == 0 && o == 4){
            #a DHC!
            if( (rdb.equiv < 1) | !("cer" %in% lois) ){
              return(NA)
            }
            n.double.bonds <-  (rdb.equiv - 1)
            n.fa.carbons <- (c - 0)
            if(n.double.bonds > max.dbl.bnds){
              return(NA)
            }
            return(paste0("Cer(h",n.fa.carbons,":",n.double.bonds,")"))
          }
          if(p == 0 && o >= 5){
            if(o == 9){
              #gal-cer
                if((rdb.equiv < 2) | !("galcer" %in% lois)){
                  return(NA)
                }
                n.double.bonds <-  (rdb.equiv - 2)
                n.fa.carbons <- (c - 6)
                
                if(n.double.bonds > max.dbl.bnds){
                  return(NA)
                }
                return(paste0("GalCer(",n.fa.carbons,":",n.double.bonds,")"))
            }
            if(o == 6){
              #S-lipid
              if((rdb.equiv < 3) | !("slipid" %in% lois)){
                return(NA)
              }
              n.double.bonds <-  (rdb.equiv - 3)
              n.fa.carbons <- (c - 3)
              
              if(n.double.bonds > max.dbl.bnds){
                return(NA)
              }
              return(paste0("S(h",n.fa.carbons,":",n.double.bonds,")"))
            }
            if(o == 5){
              #G-lipid
              if((rdb.equiv < 3) | !("glipid" %in% lois)){
                return(NA)
              }
              n.double.bonds <-  (rdb.equiv - 3)
              n.fa.carbons <- (c - 2)
              
              if(n.double.bonds > max.dbl.bnds){
                return(NA)
              }
              return(paste0("G(h",n.fa.carbons,":",n.double.bonds,")"))
            }
        
          }
          if(p == 1){
            if(o >= 9 && o <= 10 && ("gpc" %in% lois) ){
              #an GPC!
              #make sure to keep this BEFORE the PS statements, and do not put return(NA) statements here
              if((rdb.equiv >= 1 && o == 9) | (rdb.equiv >= 2 && o == 10)){
                if(o == 9){
                  #regular old GPC
                  n.double.bonds <-  (rdb.equiv - 1)
                  n.fa.carbons <- (c - 3)
                  if(n.double.bonds < max.dbl.bnds){
                    return(paste0("GPC(h",n.fa.carbons,":",n.double.bonds,")"))
                  }
                }else{  # o == 10 -> an acyl-GPC
                  n.double.bonds <-  (rdb.equiv - 2)
                  n.fa.carbons <- (c - 3)
                  if(n.double.bonds < max.dbl.bnds){
                    return(paste0("acyl-GPC(h",n.fa.carbons,":",n.double.bonds,")"))
                  }
                }

              }
            }
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
            if((o == 7 | o == 8) && ("cer1p" %in% lois)){
              #a Cer-1-P
              #make sure to keep this BEFORE the PC/PE statements, and do not put unnecessary return(NA) statements here
              if(rdb.equiv >= 1 & rdb.equiv < 2 & o == 7){
                n.double.bonds <-  (rdb.equiv - 1)
                n.fa.carbons <- (c - 0)
                if(n.double.bonds < max.dbl.bnds){
                  return(paste0("Cer-1-P(h",n.fa.carbons,":",n.double.bonds,")"))
                }
              }
              
              if(rdb.equiv >= 2 & rdb.equiv < 3 & o == 8 & strain == "p.gingivalis"){
                n.double.bonds <-  (rdb.equiv - 2)
                n.fa.carbons <- (c - 0)
                if(n.double.bonds < max.dbl.bnds){
                  return(paste0("acyl-Cer-1-P(h",n.fa.carbons,":",n.double.bonds,")"))
                }
              }
              
            }
            if(o>=7 && o<=8 ){
                #PE has 2 rdb , base 5 carbons
                if((rdb.equiv < 2) | !("pe" %in% lois)){
                #if((o == 8 && rdb.equiv < 2) | (o == 7 && rdb.equiv < 1) | !("pe" %in% lois)){
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
            if(o == 12){
              #an IPC!
              if( (rdb.equiv < 2) | !("ipc" %in% lois) ){
                return(NA)
              }
              n.double.bonds <-  (rdb.equiv - 2)
              n.fa.carbons <- (c - 6)
              if(n.double.bonds > max.dbl.bnds){
                return(NA)
              }
              return(paste0("IPC(h",n.fa.carbons,":",n.double.bonds,")"))
            }
          }
          if(length(which_to_look_for)==1){
            return(NA)
          }
        }
        if(n == 3 && s == 0){
          if(p == 1 && o == 13){
            #a diacyl-GS-PDHC!
            if( (rdb.equiv < 5) | !("dgp" %in% lois) ){
              return(NA)
            }
            n.double.bonds <-  (rdb.equiv - 5)
            n.fa.carbons <- (c - 5)
            if(n.double.bonds > max.dbl.bnds){
              return(NA)
            }
            return(paste0("GS-PDHC(h2,",n.fa.carbons,":",n.double.bonds,")"))
          }
          if(p == 1 && o == 12){
            #a acyl-GS-PDHC!
            if( (rdb.equiv < 4) | !("agp" %in% lois) ){
              return(NA)
            }
            n.double.bonds <-  (rdb.equiv - 4)
            n.fa.carbons <- (c - 5)
            if(n.double.bonds > max.dbl.bnds){
              return(NA)
            }
            return(paste0("GS-PDHC(h2,",n.fa.carbons,":",n.double.bonds,")"))
          }
          return(NA)
        }
      }
    }
  }

  return(NA)
}
most<-function(x){
  return( (length(which(x))/length(x)) > 0.95 )
}

assign_comp_from_mz <- function(mz,elements,isotope,atomic.mass,mins.list,maxs.list,rdb.rule,rdbrange,ionmode,all.elements,error.ppm,adducts){
  mz <- as.numeric(mz)
  
  if("m.minus.2h" %in% adducts){
    mz <- (2*mz)
    double.charge.multiplier <- 2
  }else{
    double.charge.multiplier <- 1
  }
  
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
      e.bound = ( 1 * double.charge.multiplier )
    }
    if(ionmode == "pos.ion"){
      e.bound = -1
    }

    bounds <- list(lower = list(ind = c(1:length(v)), val = c(mins.list,e.bound)),
                   upper = list(ind = c(1:length(v)), val = c(maxs.list,e.bound)))
    
    rdb.row <- rep(0,length(v))
    rdb.row[ele.list %in% c("H1","D2","T3","Cl35","Cl37","Na23","Li7","Li6")] <- (-0.5)
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
                       (a$optimum/double.charge.multiplier),
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
                       (b$optimum/double.charge.multiplier),
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
      e.bound = ( 1 * double.charge.multiplier )
      
    }
    if(ionmode == "pos.ion"){
      e.bound = -1
    }
    
    min.r = floor((rhs.var-rdbrange[2])/2)-1
    max.r = ceiling((rhs.var-rdbrange[1])/2)+1

    bounds <- list(lower = list(ind = c(1:length(v)), val = c(mins.list,e.bound,min.r)),
                   upper = list(ind = c(1:length(v)), val = c(maxs.list,e.bound,max.r)))
    
    rdb.row <- rep(0,length(v))
    rdb.row[ele.list %in% c("H1","D2","T3","Cl35","Cl37","Na23","Li7","Li6")] <- (-0.5)
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
                 (a$optimum/double.charge.multiplier),
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
                 (b$optimum/double.charge.multiplier),
                 ((b$optimum-x)/x)*10^6))
      }
    }else{
      return(c(rep("",3)))
    }}))
    
  }
  
}

calculate_logical_element_minimums <- function(loi,lipidclasslist,elements,isotopes,current.mins){

  mins <- current.mins 
  
  mins[paste0(elements,isotopes) == "C12"] <- min(atomic_min_max_table$min.c.12[lipidclasslist %in% loi]) 
  mins[paste0(elements,isotopes) == "H1"] <- min(atomic_min_max_table$min.h.1[lipidclasslist %in% loi])
  mins[paste0(elements,isotopes) == "O16"] <- min(atomic_min_max_table$min.o.16[lipidclasslist %in% loi]) 
  mins[paste0(elements,isotopes) == "N14"] <- min(atomic_min_max_table$min.n.14[lipidclasslist %in% loi])
  mins[paste0(elements,isotopes) == "P31"] <- min(atomic_min_max_table$min.p.31[lipidclasslist %in% loi]) 
  mins[paste0(elements,isotopes) == "S32"] <- min(atomic_min_max_table$min.s.32[lipidclasslist %in% loi]) 
  
  return(mins)
  
}
calculate_logical_element_maximums <- function(loi,lipidclasslist,elements,isotopes,current.maxs){

  maxs <- current.maxs 
  
  maxs[paste0(elements,isotopes) == "C12"] <- max(atomic_min_max_table$max.c.12[lipidclasslist %in% loi]) 
  maxs[paste0(elements,isotopes) == "H1"] <- max(atomic_min_max_table$max.h.1[lipidclasslist %in% loi]) 
  maxs[paste0(elements,isotopes) == "O16"] <- max(atomic_min_max_table$max.o.16[lipidclasslist %in% loi])
  maxs[paste0(elements,isotopes) == "N14"] <- max(atomic_min_max_table$max.n.14[lipidclasslist %in% loi])
  maxs[paste0(elements,isotopes) == "P31"] <- max(atomic_min_max_table$max.p.31[lipidclasslist %in% loi])
  maxs[paste0(elements,isotopes) == "S32"] <- max(atomic_min_max_table$max.s.32[lipidclasslist %in% loi])
  
  return(maxs)
  
}

# make lists
lipidclasslist <- c("DAG" = "dag","TAG" = "tag","PA" = "pa","PG" = "pg","PI" = "pi","CL" = "cl","SM" = "sm","Ceramide" = "cer","PC" = "pc","PE" = "pe","PS" = "ps","FFA" = "ffa","GDG" = "gdg",
                    "IPC" = "ipc","SPC" = "spc","EPC" = "epc", "GPC" = "gpc", "Cer-1-P" = "cer1p","Gal-Cer" = "galcer","S-lipid" = "slipid",
                    "G-lipid" = "glipid","GS-lipid" = "gslipid","GS-PA" = "gspa","DHC-PIP-DHC" = "dpd","diacyl-GS-PDHC" = "dgp","acyl-GS-PDHC" = "agp","GS-PG" = "gspg","Sulfatides" = "sulf","Gangliosides" = "gng")
negionclasses <- c("pe","cer","cl","pi","pg","pa","ps","ffa","gdg","ipc","spc","epc","cer1p","gpc","galcer","slipid","glipid","gslipid","gspa","dpd", #"dgp","agp",
                   "gspg","sulf","gng")
posionclasses <- c("pc","sm","tag","dag","gdg")
negion_lipids_in_both_baceriodetes_and_eukaryotes <- c("pe","ps","galcer")
bacteroidetes_lipids <- c("ipc","spc","epc","gpc","cer1p","slipid","glipid","gslipid","cer","gspa","dpd","pe","ps","gspg")
lipidclasslist_fullnames <- c("Cardiolipin" = "CL", "Ceramide-1-phosphate" = "Cer-1-P", "Diacylglycerol" = "DAG", "Dihydroceramide(DHC)-phosphoinositolphosphoryl-DHC" = "DHC-PIP-DHC", "Ethanolamine phosphoryl ceramide" = "EPC", "Free fatty acids" = "FFA", "Galactosylceramide" = "Gal-Cer","Glycerol phosphoryl ceramide" = "GPC", "Glycine lipid" = "G-lipid", "Glycosyldiacylglycerols" = "GDG", "Glycylserine-phosphoglycerol" = "GS-PG", "Glycylserine-phosphoryldihydroceramide" = "GS-PDHC",
                              "Inositol phosphorylceramide" = "IPC", "Lipoglycylserine" = "GS-lipid", "Lipoglycylserine phosphatidic acid" = "GS-PA", "Phosphatidic acid" = "PA", "Phosphatidylcholine" = "PC", "Phosphatidylethanolamine" = "PE","Phosphatidylglycerol" = "PG",  "Phosphatidylinositol" = "PI","Phosphatidylserine" = "PS",  "Serine lipids" = "S-lipid","Serine phosphoryl ceramide" = "SPC",  "Sphingomyelin" = "SM",   "Triacylglycerol" = "TAG")
lipidclasslist_fullnames <- lipidclasslist_fullnames[order(lipidclasslist_fullnames)]

relevant.col.names <- c("mz","intensity","rel.intensity","theo.mass","delta","rdb.equiv","composition")

pls <- read.csv('pl_database.csv')

strainslist <- c("n/a",unique(pls$strain))
strainslist <- strainslist[strainslist != ""]
names(strainslist) <- strainslist

specieslist <- c("n/a",unique(pls$species))
specieslist <- specieslist[specieslist != ""]
names(specieslist) <- specieslist

atomic_min_max_table <- read.csv('atomic_min_max_values.csv')
atomic_min_max_table <- atomic_min_max_table[match(lipidclasslist,atomic_min_max_table$species),]

all.elements = read.csv('atomics_wts_and_isotopic_compositions_for_all_elements.csv')
all.elements$full.element = paste0(all.elements$element,"[",all.elements$isotope,"]")
this_table = all.elements[paste0(all.elements$element,all.elements$isotope) %in% c("H1","C12","N14","P31","O16","S32"),
                          colnames(all.elements) %in% c("element","isotope","atomic.mass")]
this_table$minimum <- rep(0,nrow(this_table))
this_table$maximum <- rep(150,nrow(this_table))
rownames(this_table) <- NULL
this_table2 <- this_table
choices = paste0(all.elements$element,"[",all.elements$isotope,"]")
names(choices) = all.elements$full.element

# these make the css and js for the tooltip dropdown that shows all of the lipid classes and their full names
css <- '
.tooltip {
  pointer-events: none;
}
.tooltip > .tooltip-inner {
  pointer-events: none;
  background-color: rgba(64, 64, 64, 0.95); /* Dark gray, partially opaque background */
  color: #FFFFFF; /* White text */
  border: 1px solid #FFFFFF; /* White border */
  padding: 10px;
  font-size: 14px; /* Adjusted font size */
  text-align: left; /* Align text to the left */
  margin-left: 0;
  max-width: 1000px; /* Remove the maximum width */
  width: 500px; /* Set a specific width if needed */
  font-weight: normal; /* Set to normal weight (not bold) */
  font-style: normal; /* Set to normal style (not italicized) */
  html: TRUE;
}
.tooltip > .arrow::before {
  border-right-color: rgba(64, 64, 64, 0.8); /* Dark gray, partially opaque background for arrow */
}
'

js <- "
$(function () {
  $('[data-toggle=tooltip]').tooltip()
})
"

ui <- fluidPage(
  tags$head(
    tags$style(HTML(css)),
    tags$script(HTML(js))
  ),
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
                   condition = "input.tabsid == 0  && input.domain == 'euk'",
                   selectInput("species","Eukaryotic species:",
                               specieslist,
                               selected = "n/a",
                               multiple = F)
                 ),
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
                                        "[M+Na]+" = "m.plus.sodium",
                                        "[M+Li]+" = "m.plus.lithium")),
                   selectInput("posionloi",span("Lipid classes to search for:",
                                                # https://stackoverflow.com/questions/73306444/r-shiny-popup-window-when-hovering-over-icon
                                                # https://stackoverflow.com/questions/75515925/r-shiny-multiple-onhover-popups-each-with-different-css-style
                                                # https://stackoverflow.com/questions/3340802/add-line-break-within-tooltips
                                                div(`data-toggle` = "tooltip", 
                                                    `data-placement` = "right",
                                                    `data-html`="true",
                                                    style = "display:inline-block;",
                                                    title = paste(mapply(function(x,y) paste0(x," = ",y),lipidclasslist_fullnames,names(lipidclasslist_fullnames)),collapse = "<br/>"),
                                                    icon("info-circle"))
                   ),
                               lipidclasslist,
                               selected = posionclasses[!(posionclasses %in% bacteroidetes_lipids)],
                               multiple = T)
                 ),
                 conditionalPanel(
                   condition = "input.tabsid == 0  && input.ionmode == 'neg.ion'",
                   checkboxGroupInput("negionadduct", "Ions/adducts to look for:",
                                      c("[M-H]-" = "m.minus.h",
                                        "[M+Cl]-" = "m.plus.chloride",
                                        "[M-2H]2-" = "m.minus.2h")),
                   selectInput("negionloi", span("Lipid classes to search for:",
                                                  div(`data-toggle` = "tooltip", 
                                                      `data-placement` = "right",
                                                      `data-html`="true",
                                                      style = "display:inline-block;",
                                                      title = paste(mapply(function(x,y) paste0(x," = ",y),lipidclasslist_fullnames,names(lipidclasslist_fullnames)),collapse = "<br/>"),
                                                      icon("info-circle"))
                                                  ),
                                 lipidclasslist,
                                 selected = c( (negionclasses[!(negionclasses %in% bacteroidetes_lipids)]) ,negion_lipids_in_both_baceriodetes_and_eukaryotes),
                                 multiple = T),
                   
                 ),
                 conditionalPanel(
                   condition = "input.tabsid == 0  && input.ionmode == 'neutral'",
                   selectInput("neutralloi",span("Lipid classes to search for:",
                                                 div(`data-toggle` = "tooltip", 
                                                     `data-placement` = "right",
                                                     `data-html`="true",
                                                     style = "display:inline-block;",
                                                     title = paste(mapply(function(x,y) paste0(x," = ",y),lipidclasslist_fullnames,names(lipidclasslist_fullnames)),collapse = "<br/>"),
                                                     icon("info-circle"))
                   ),
                               lipidclasslist,
                               selected = c( (c(negionclasses,posionclasses)[!(c(negionclasses,posionclasses) %in% bacteroidetes_lipids)]) , negion_lipids_in_both_baceriodetes_and_eukaryotes),
                               multiple = T)
                 ),
                 sliderInput("maxdblbnds",
                             "Maximum number of acyl/alkyl double bonds allowed",
                             min = 0,max = 50,value = 12,step = 1,round = T,animate = F, ticks = F),
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
                                                min = -5, 
                                                max = 30, 
                                                value = c(-1,25),step = 0.5 ,ticks = F  )),
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
                        actionButton("add_btn", "Add/Update element",icon("check"),style="color: #FFFFFF; background-color: #7DA8F2; border-color: #7DA8F2"),
                        actionButton("delete_btn", "Delete element",icon("trash"),style="color: #FFFFFF; background-color: #FF0000; border-color: #2e6da4"),
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
                   condition = "input.tabsid == 1  && input.bulkdomain == 'euk'",
                   selectInput("bulkspecies","Eukaryotic species:",
                               specieslist,
                               selected = "n/a",
                               multiple = F)
                 ),
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
                                        "[M+Na]+" = "m.plus.sodium",
                                        "[M+Li]+" = "m.plus.lithium")),
                   selectInput("bulkposionloi",span("Lipid classes to search for:",
                                                    div(`data-toggle` = "tooltip", 
                                                        `data-placement` = "right",
                                                        `data-html`="true",
                                                        style = "display:inline-block;",
                                                        title = paste(mapply(function(x,y) paste0(x," = ",y),lipidclasslist_fullnames,names(lipidclasslist_fullnames)),collapse = "<br/>"),
                                                        icon("info-circle"))
                   ),
                               lipidclasslist,
                               selected = posionclasses[!(posionclasses %in% bacteroidetes_lipids)],
                               multiple = T)
                 ),
                 conditionalPanel(
                   condition = "input.tabsid == 1  && input.bulkionmode == 'neg.ion'",
                   checkboxGroupInput("bulknegionadduct", "Ions/adducts to look for:",
                                      c("[M-H]-" = "m.minus.h",
                                        "[M+Cl]-" = "m.plus.chloride",
                                        "[M-2H]2-" = "m.minus.2h")),
                   selectInput("bulknegionloi",span("Lipid classes to search for:",
                                                    div(`data-toggle` = "tooltip", 
                                                        `data-placement` = "right",
                                                        `data-html`="true",
                                                        style = "display:inline-block;",
                                                        title = paste(mapply(function(x,y) paste0(x," = ",y),lipidclasslist_fullnames,names(lipidclasslist_fullnames)),collapse = "<br/>"),
                                                        icon("info-circle"))
                   ),
                               lipidclasslist,
                               selected = c( (negionclasses[!(negionclasses %in% bacteroidetes_lipids)]) ,negion_lipids_in_both_baceriodetes_and_eukaryotes),
                               multiple = T)
                 ),
                 conditionalPanel(
                   condition = "input.tabsid == 1  && input.bulkionmode == 'neutral'",
                   selectInput("bulkneutralloi",span("Lipid classes to search for:",
                                                     div(`data-toggle` = "tooltip", 
                                                         `data-placement` = "right",
                                                         `data-html`="true",
                                                         style = "display:inline-block;",
                                                         title = paste(mapply(function(x,y) paste0(x," = ",y),lipidclasslist_fullnames,names(lipidclasslist_fullnames)),collapse = "<br/>"),
                                                         icon("info-circle"))
                   ),
                               lipidclasslist,
                               selected = c( (c(negionclasses,posionclasses)[!(c(negionclasses,posionclasses) %in% bacteroidetes_lipids) ]) ,negion_lipids_in_both_baceriodetes_and_eukaryotes),
                               multiple = T)
                 ),
                 radioButtons("returnAll", "Would you like to...", 
                              choiceNames = c("Return all peaks",
                                              "Only return peaks with a species assigned"),
                              choiceValues = c("returnallpeaks","returnassignedpeaks"),
                              inline = F, selected = ("returnallpeaks")),
                 sliderInput("bulkmaxdblbnds",
                             "Max number of acyl/alkyl double bonds allowed",
                             min = 0,max = 50,value = 12,step = 1,round = T,animate = F, ticks = F),
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
                        h1(strong("4. Input data:"),downloadLink("testdataset1", label = "(example data)") , style = "font-size:22.5px;"),
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
                                                               min = -5, 
                                                               max = 30, 
                                                               value = c(-1,25),step = 0.5,ticks = F)),
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
                                       actionButton("bulk.add_btn", "Add/Update element",icon("check"),style="color: #FFFFFF; background-color: #7DA8F2; border-color: #7DA8F2"),
                                       actionButton("bulk.delete_btn", "Delete element",icon("trash"),style="color: #FFFFFF; background-color: #FF0000; border-color: #2e6da4"),
                                       DTOutput("bulk_shiny_table",width = "20%"),
                                       conditionalPanel(condition = "input.tabsid == 1 && input.transferdata == 'y' ",
                                                        h1(strong("4. Recalibrated data has been automatically imported from the recalibration tab"), style = "font-size:22.5px;"),
                                                        radioButtons("transferdata_undo",
                                                                     "Would you like to proceed with the data from the recalibration tab, or would you like to enter new data in this tab instead?",
                                                                     choices = c("Proceed with recalibrated data" = "y",
                                                                                 "Input new data in this tab" = "n"),
                                                                     selected = "y",
                                                                     inline = T)),
                                       conditionalPanel(condition = "input.tabsid == 1 && input.transferdata !== 'y' ",
                                                        h1(strong("4. Input data:"),downloadLink("testdataset2", label = "(example data)"),style = "font-size:22.5px;"),
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
                                    actionButton("bulk.analyze.button", "Analyze data", icon("redo"),style="color: #FFFFFF; background-color: #5AB26D; border-color: #5AB26D; margin-top:25px",width = '100%')
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
                 h1(strong("1. Input data:"),downloadLink("testdataset3", label = "(example data)"), style = "font-size:22.5px;"),
                 import_copypaste_ui3("myid3",
                                      title = ""),
                 tags$b("Import status:"),
                 htmlOutput(outputId = "status3"),
                 tags$b("Sample name:"),
                 verbatimTextOutput(outputId = "name3")
               ),
               column(width = 8,
                      h1(strong("2. Select peak(s) you wish to use to recalibrate:"), style = "font-size:22.5px;"),
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
                                  downloadButton("download3", "Download aligned data",style="color: #FFFFFF; background-color: #7DA8F2; border-color: #7DA8F2")
                                  )
                      )
             )),
    )
  
)

server <- function(input, output, session){
  #add chlorine to the elemental table when user adds chlorine adduct
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
    if("m.minus.2h" %in% input$negionadduct){
      f.e <- "C[13]"
      t <- data.frame(element = all.elements$element[match(f.e,all.elements$full.element)],
                      isotope = all.elements$isotope[match(f.e,all.elements$full.element)],
                      atomic.mass = all.elements$atomic.mass[match(f.e,all.elements$full.element)],
                      minimum = 0,
                      maximum = 1)
      if( !(paste0(t$element,t$isotope) %in% paste0(this_table()$element,this_table()$isotope)) ){
        this_table(rbind(t, this_table()))
      }
    }else{
      element_is_carbon <- (this_table()$element == "C")
      isotope_is_thirteen <- (this_table()$isotope == "13")
      this_table(this_table()[!(element_is_carbon & isotope_is_thirteen),])
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
    if("m.minus.2h" %in% input$bulknegionadduct){
      f.e <- "C[13]"
      t <- data.frame(element = all.elements$element[match(f.e,all.elements$full.element)],
                      isotope = all.elements$isotope[match(f.e,all.elements$full.element)],
                      atomic.mass = all.elements$atomic.mass[match(f.e,all.elements$full.element)],
                      minimum = 0,
                      maximum = 1)
      if( !(paste0(t$element,t$isotope) %in% paste0(this_table2()$element,this_table2()$isotope)) ){
        this_table2(rbind(t, this_table2()))
      }
    }else{
      element_is_carbon <- (this_table2()$element == "C")
      isotope_is_thirteen <- (this_table2()$isotope == "13")
      this_table2(this_table2()[!(element_is_carbon & isotope_is_thirteen),])
    }
  })
  
  #add sodium to the elemental table when user adds sodium adduct
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
    
    if("m.plus.lithium" %in% input$posionadduct){
      f.e <- "Li[7]"
      t <- data.frame(element = all.elements$element[match(f.e,all.elements$full.element)],
                      isotope = all.elements$isotope[match(f.e,all.elements$full.element)],
                      atomic.mass = all.elements$atomic.mass[match(f.e,all.elements$full.element)],
                      minimum = 0,
                      maximum = 1)
      if( !(paste0(t$element,t$isotope) %in% paste0(this_table()$element,this_table()$isotope)) ){
        this_table(rbind(t, this_table()))
      }
    }else{
      this_table(this_table()[this_table()$element != "Li",])
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
    
    if("m.plus.lithium" %in% input$bulkposionadduct){
      f.e <- "Li[7]"
      t <- data.frame(element = all.elements$element[match(f.e,all.elements$full.element)],
                      isotope = all.elements$isotope[match(f.e,all.elements$full.element)],
                      atomic.mass = all.elements$atomic.mass[match(f.e,all.elements$full.element)],
                      minimum = 0,
                      maximum = 1)
      if( !(paste0(t$element,t$isotope) %in% paste0(this_table2()$element,this_table2()$isotope)) ){
        this_table2(rbind(t, this_table2()))
      }
    }else{
      this_table2(this_table2()[this_table2()$element != "Li",])
    }
  })
  
  #make it so that selecting [M-2H]2- de-selects other ions and vice versa
  temporary_lipid_negionadduct_list <- reactiveVal()
  observeEvent(input$negionadduct,{
    if(length(input$negionadduct) > 1){
      if("m.minus.2h" %in% input$negionadduct){ #more than one lipid selected, and one of them is m.minus.2h -> we have something to worry about
        #first, identify the newly selected adduct by finding the adduct in input$negionadduct which is not in temporary_lipid_negionadduct_list()
        new_lipid <- input$negionadduct[!(input$negionadduct %in% temporary_lipid_negionadduct_list())]
        updateSelectInput(session = session,
                          "negionadduct",
                          selected = new_lipid)
        temporary_lipid_negionadduct_list(new_lipid)
      }else{
        #more than one adduct selected, but none of them are m.minus.2h -> nothing to worry about
        temporary_lipid_negionadduct_list(input$negionadduct)
      }
    }else{
      # 0-1 adduct is selected -> nothing to worry about
      temporary_lipid_negionadduct_list(input$negionadduct)
    }
  })
  temporary_lipid_bulknegionadduct_list <- reactiveVal()
  observeEvent(input$bulknegionadduct,{
    if(length(input$bulknegionadduct) > 1){
      if("m.minus.2h" %in% input$bulknegionadduct){ #more than one lipid selected, and one of them is m.minus.2h -> we have something to worry about
        #first, identify the newly selected adduct by finding the adduct in input$bulknegionadduct which is not in temporary_lipid_bulknegionadduct_list()
        new_lipid <- input$bulknegionadduct[!(input$bulknegionadduct %in% temporary_lipid_bulknegionadduct_list())]
        updateSelectInput(session = session,
                          "bulknegionadduct",
                          selected = new_lipid)
        temporary_lipid_bulknegionadduct_list(new_lipid)
      }else{
        #more than one adduct selected, but none of them are m.minus.2h -> nothing to worry about
        temporary_lipid_bulknegionadduct_list(input$bulknegionadduct)
      }
    }else{
      # 0-1 adduct is selected -> nothing to worry about
      temporary_lipid_bulknegionadduct_list(input$bulknegionadduct)
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
  
  #reset species to n/a everytime domain is changed
  observeEvent(input$domain,
               updateSelectInput(session = session,
                                 "species",
                                 "Eukaryotic species:",
                                 choices = specieslist,
                                 selected = "n/a"
               )
  )
  observeEvent(input$bulkdomain,
               updateSelectInput(session = session,
                                 "bulkspecies",
                                 "Eukaryotic species:",
                                 choices = specieslist,
                                 selected = "n/a"
               )
  )
  
  #make it so that when user says "no" in the second tab to not transfer data from the recal tab, that it makes input$transferdata == "no", and vice versa
  observeEvent(input$transferdata_undo,{
    if(input$transferdata_undo == "n"){
      updateRadioButtons(session = session,
                         inputId = "transferdata",
                         "Transfer recalibrated data to second tab?",
                         choices = c("Yes" = "y",
                                     "No" = "n"),
                         selected = "n",
                         inline = T
      )
      
    }
  }
  )
  observeEvent(input$transferdata,{
    if(input$transferdata == "y"){
      updateRadioButtons(session = session,
                         inputId = "transferdata_undo",
                         "Would you like to proceed with the data from the recalibration tab, or would you like to enter new data in this tab instead?",
                         choices = c("Proceed with recalibrated data" = "y",
                                     "Input new data in this tab" = "n"),
                         selected = "y",
                         inline = T
      )
      
    }
  }
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
          
                 #clear all negative ion mode ions/adducts 
                 updateCheckboxGroupInput(session = session,
                                          "negionadduct",
                                          "Ions/adducts to look for:",
                                          c("[M-H]-" = "m.minus.h",
                                            "[M+Cl]-" = "m.plus.chloride",
                                            "[M-2H]2-" = "m.minus.2h"),
                                          selected=character(0))
                 if(input$strain %in% c("bacteroidetes","p.gingivalis") ){
                   updateSelectInput(session = session, "negionloi",
                                     selected = negionclasses[(negionclasses %in% bacteroidetes_lipids)])
                   updateSelectInput(session = session, "neutralloi",
                                     selected = c(negionclasses,posionclasses)[(c(negionclasses,posionclasses) %in% bacteroidetes_lipids)])
                 }else{
                   updateSelectInput(session = session, "negionloi",
                                     selected = c( negionclasses[!(negionclasses %in% bacteroidetes_lipids)] , negion_lipids_in_both_baceriodetes_and_eukaryotes) )
                   updateSelectInput(session = session, "neutralloi",
                                     selected = c(c(negionclasses,posionclasses)[!(c(negionclasses,posionclasses) %in% bacteroidetes_lipids) ],negion_lipids_in_both_baceriodetes_and_eukaryotes) )
                 }
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
                                              "[M+Na]+" = "m.plus.sodium",
                                              "[M+Li]+" = "m.plus.lithium"),
                                            selected=character(0))
                   if(input$strain %in% c("bacteroidetes","p.gingivalis") ){
                     updateSelectInput(session = session, "posionloi",
                                       selected = posionclasses)
                     updateSelectInput(session = session, "neutralloi",
                                       selected = c(negionclasses,posionclasses)[(c(negionclasses,posionclasses) %in% bacteroidetes_lipids)])
                   }else{
                     updateSelectInput(session = session, "posionloi",
                                       selected = posionclasses[!(posionclasses %in% bacteroidetes_lipids)])
                     updateSelectInput(session = session, "neutralloi",
                                       selected = c(c(negionclasses,posionclasses)[!(c(negionclasses,posionclasses) %in% bacteroidetes_lipids) ],negion_lipids_in_both_baceriodetes_and_eukaryotes) )
                   }

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
                   updateCheckboxGroupInput(session = session,
                                            "posionadduct",
                                            "Ions/adducts to look for:",
                                            c("[M+H]+" = "m.plus.h",
                                              "[M+NH4]+" = "m.plus.ammonia",
                                              "[M+Na]+" = "m.plus.sodium",
                                              "[M+Li]+" = "m.plus.lithium"),
                                            selected=character(0))
                   if(input$strain %in% c("bacteroidetes","p.gingivalis") ){
                     updateSelectInput(session = session, "negionloi",
                                       selected = negionclasses[(negionclasses %in% bacteroidetes_lipids)])
                     updateSelectInput(session = session, "posionloi",
                                       selected = posionclasses)
                   }else{
                     updateSelectInput(session = session, "negionloi",
                                       selected = c(negionclasses[!(negionclasses %in% bacteroidetes_lipids)],negion_lipids_in_both_baceriodetes_and_eukaryotes))
                     updateSelectInput(session = session, "posionloi",
                                       selected = posionclasses[!(posionclasses %in% bacteroidetes_lipids)])
                   }

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
                 if(input$bulkstrain %in% c("bacteroidetes","p.gingivalis") ){
                   updateSelectInput(session = session, "bulknegionloi",
                                     selected = negionclasses[(negionclasses %in% bacteroidetes_lipids)])
                   updateSelectInput(session = session, "bulkneutralloi",
                                     selected = c(negionclasses,posionclasses)[(c(negionclasses,posionclasses) %in% bacteroidetes_lipids)])
                 }else{
                   updateSelectInput(session = session, "bulknegionloi",
                                     selected = c(negionclasses[!(negionclasses %in% bacteroidetes_lipids)],negion_lipids_in_both_baceriodetes_and_eukaryotes) )
                   updateSelectInput(session = session, "bulkneutralloi",
                                     selected = c(c(negionclasses,posionclasses)[!(c(negionclasses,posionclasses) %in% bacteroidetes_lipids) ],negion_lipids_in_both_baceriodetes_and_eukaryotes) )
                 }

                 
               }
                 if(input$bulkionmode == "neg.ion"){
                   updateCheckboxGroupInput(session = session,
                                            "bulkposionadduct",
                                            "Ions/adducts to look for:",
                                            c("[M+H]+" = "m.plus.h",
                                              "[M+NH4]+" = "m.plus.ammonia",
                                              "[M+Na]+" = "m.plus.sodium",
                                              "[M+Li]+" = "m.plus.lithium"),
                                            selected=character(0))
                   if(input$bulkstrain %in% c("bacteroidetes","p.gingivalis") ){
                     updateSelectInput(session = session, "bulkposionloi",
                                       selected = posionclasses)
                     updateSelectInput(session = session, "bulkneutralloi",
                                       selected = c(negionclasses,posionclasses)[(c(negionclasses,posionclasses) %in% bacteroidetes_lipids)])
                   }else{
                     updateSelectInput(session = session, "bulkposionloi",
                                       selected = posionclasses[!(posionclasses %in% bacteroidetes_lipids)])
                     updateSelectInput(session = session, "bulkneutralloi",
                                       selected = c(c(negionclasses,posionclasses)[!(c(negionclasses,posionclasses) %in% bacteroidetes_lipids) ],negion_lipids_in_both_baceriodetes_and_eukaryotes) )
                     
                   }
                   
                 }
                 if(input$bulkionmode == "neutral"){
                   updateCheckboxGroupInput(session = session,
                                            "bulknegionadduct",
                                            "Ions/adducts to look for:",
                                            c("[M-H]-" = "m.minus.h",
                                              "[M+Cl]-" = "m.plus.chloride",
                                              "[M-2H]2-" = "m.minus.2h"),
                                            selected=character(0))
                   updateCheckboxGroupInput(session = session,
                                            "bulkposionadduct",
                                            "Ions/adducts to look for:",
                                            c("[M+H]+" = "m.plus.h",
                                              "[M+NH4]+" = "m.plus.ammonia",
                                              "[M+Na]+" = "m.plus.sodium",
                                              "[M+Li]+" = "m.plus.lithium"),
                                            selected=character(0))
                   if(input$bulkstrain %in% c("bacteroidetes","p.gingivalis") ){
                     updateSelectInput(session = session, "bulknegionloi",
                                       selected = negionclasses[(negionclasses %in% bacteroidetes_lipids)])
                     updateSelectInput(session = session, "bulkposionloi",
                                       selected = posionclasses)
                   }else{
                     updateSelectInput(session = session, "bulknegionloi",
                                       selected = c(negionclasses[!(negionclasses %in% bacteroidetes_lipids)],negion_lipids_in_both_baceriodetes_and_eukaryotes))
                     updateSelectInput(session = session, "bulkposionloi",
                                       selected = posionclasses[!(posionclasses %in% bacteroidetes_lipids)])
                   }
                 }
               },
               ignoreInit = T)
  
  #update lipids of interest if user selects bacteriodetes strain
  observeEvent(input$strain,{
    if(input$strain %in% c("bacteroidetes","p.gingivalis") ){
      updateSelectInput(session = session, "negionloi",
                        selected = negionclasses[(negionclasses %in% bacteroidetes_lipids)])
      updateSelectInput(session = session, "posionloi",
                        selected = posionclasses)
      updateSelectInput(session = session, "neutralloi",
                        selected = c(negionclasses,posionclasses)[(c(negionclasses,posionclasses) %in% bacteroidetes_lipids)])
    }else{
      updateSelectInput(session = session, "negionloi",
                        selected = c(negionclasses[!(negionclasses %in% bacteroidetes_lipids)],negion_lipids_in_both_baceriodetes_and_eukaryotes))
      updateSelectInput(session = session, "posionloi",
                        selected = posionclasses[!(posionclasses %in% bacteroidetes_lipids)])
      updateSelectInput(session = session, "neutralloi",
                        selected = c(c(negionclasses,posionclasses)[!(c(negionclasses,posionclasses) %in% bacteroidetes_lipids) ],negion_lipids_in_both_baceriodetes_and_eukaryotes))
    }
               },ignoreInit = T)
  observeEvent(input$bulkstrain,{
    if(input$bulkstrain %in% c("bacteroidetes","p.gingivalis") ){
      updateSelectInput(session = session, "bulknegionloi",
                        selected = negionclasses[(negionclasses %in% bacteroidetes_lipids)])
      updateSelectInput(session = session, "bulkposionloi",
                        selected = posionclasses)
      updateSelectInput(session = session, "bulkneutralloi",
                        selected = c(negionclasses,posionclasses)[(c(negionclasses,posionclasses) %in% bacteroidetes_lipids)])
    }else{
      updateSelectInput(session = session, "bulknegionloi",
                        selected = c(negionclasses[!(negionclasses %in% bacteroidetes_lipids)],negion_lipids_in_both_baceriodetes_and_eukaryotes))
      updateSelectInput(session = session, "bulkposionloi",
                        selected = posionclasses[!(posionclasses %in% bacteroidetes_lipids)])
      updateSelectInput(session = session, "bulkneutralloi",
                        selected = c(c(negionclasses,posionclasses)[!(c(negionclasses,posionclasses) %in% bacteroidetes_lipids) ],negion_lipids_in_both_baceriodetes_and_eukaryotes))
    }
  },ignoreInit = T)
  
  #make it so that if user has input valid data to be realigned, tab 2 automatically goes to "from m/z" mode
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
    
    if(input$ionmode == "pos.ion"){
      adducts <- input$posionadduct
      lois <- input$posionloi
    }
    if(input$ionmode == "neg.ion"){
      adducts <- input$negionadduct
      lois <- input$negionloi
    }
    if(input$ionmode == "neutral"){
      
      # I removed m.minus.2h from adducts list for neutral ion mode on 11/17/23
      # It makes the least sense and would take a lot of work to add...
      
      adducts <- c("m.plus.h","m.plus.ammonia","m.plus.sodium","m.minus.h","m.plus.chloride","m.plus.lithium")
      lois <- input$neutralloi
    }
    
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
                                         error.ppm = input$delta.error,
                                         adducts = adducts)
      comp <- unlist(lapply(comp.object,`[[`, 1))
      theo.mass <- as.numeric(unlist(lapply(comp.object,`[[`, 2)))
      
      #extract numbers of each element
      c <- extract_num_elements("C",comp) + extract_num_elements("C[13]",comp) + extract_num_elements("C[14]",comp)
      h <- extract_num_elements("H",comp) + extract_num_elements("D[2]",comp) + extract_num_elements("T[3]",comp)
      o <- extract_num_elements("O",comp) + extract_num_elements("O[17]",comp) + extract_num_elements("O[18]",comp)
      n <- extract_num_elements("N",comp) + extract_num_elements("N[15]",comp)
      p <- extract_num_elements("P",comp)
      s <- extract_num_elements("S",comp) + extract_num_elements("S[33]",comp) + extract_num_elements("S[34]",comp) + extract_num_elements("S[36]",comp)
      na <- extract_num_elements("Na",comp)
      cl <- extract_num_elements("Cl[35]",comp) + extract_num_elements("Cl[37]",comp)
      li <- extract_num_elements("Li",comp) + extract_num_elements("Li[6]",comp)
    }else{
      #input is not numeric
      c <- extract_num_elements("C",input$comp)
      h <- extract_num_elements("H",input$comp)
      o <- extract_num_elements("O",input$comp)
      n <- extract_num_elements("N",input$comp)
      p <- extract_num_elements("P",input$comp)
      s <- extract_num_elements("S",input$comp) + extract_num_elements("S[33]",input$comp) + extract_num_elements("S[34]",input$comp) + extract_num_elements("S[36]",input$comp)
      na <- extract_num_elements("Na",input$comp)
      cl <- extract_num_elements("Cl",input$comp)
      li <- extract_num_elements("Li",input$comp) + extract_num_elements("Li[6]",input$comp)
    }
    
    str1 <- paste0("Carbon: ",c)
    str2 <- paste0("Hydrogen: ",h)
    str3 <- paste0("Oxygen: ",o)
    str4 <- paste0("Nitrogen: ",n)
    str5 <- paste0("Phosphorus: ",p)
    str5_5 <- paste0("Sulfur: ",s)
    str6 <- paste0("Sodium: ",na)
    str7 <- paste0("Chlorine: ",cl)
    str7_5 <- paste0("Lithium: ",li)
    str8<-""
    str9<-""
    rdb.equiv = ( c - ((h+cl+na+li)/2) + ((n+p)/2) +1 )
    str10<-paste0("RDB equiv. = ", rdb.equiv)
    str11<-""
    
    species.result <- assign_species(c,h,o,n,p,s,na,cl,li,input$ionmode,adducts,rdb.equiv,input$domain,lois,input$maxdblbnds,input$strain)
    if(!is.na(species.result)){
      if(species.result %in% pls$gen.structure){
        
        if( (input$strain != "n/a" && species.result %in% pls[pls$strain == input$strain , colnames(pls) == "gen.structure"]) | (input$species != "n/a" && species.result %in% pls[pls$species == input$species , colnames(pls) == "gen.structure"]) ){
          
          if(input$strain != "n/a" && species.result %in% pls[pls$strain == input$strain , colnames(pls) == "gen.structure"]){
            pls2 <- pls[pls$strain == input$strain,]
          }else{ #then the following bool must return T: 'input$species != "n/a" && species.result %in% pls[pls$species == input$species , colnames(pls) == "gen.structure"]'
            pls2 <- pls[pls$species == input$species,]
          }
          
          c2 <- unname(unlist(pls2[ match(species.result,pls2$gen.structure) , 2:ncol(pls2) ]))
          c2 <- c2[c2 != ""]
          c2 <- c2[-length(c2)]
          if(input$strain != "n/a"){
            c2 <- c(c2,paste0(" (",input$strain," specific)"))
          }else{
            c2 <- c(c2,paste0(" (",input$species," specific)"))
          }
          
          if(length(c2) <= 2){
            species.result <- paste0(species.result ,"</b>","; most likely exact structure: ", paste0(c2,collapse = "") )
          }else{
            species.result <- paste0(species.result ,"</b>","; most likely exact structures: ", paste(c2, collapse = ", ") )
          }
        }else{
          pls2 <- pls[pls$strain == "" & pls$species == "",]
          
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
      str5_5 <- ""
      str6 <- ""
      str7 <- paste0("Delta is: ",format(round(((theo.mass-as.numeric(input$comp))/theo.mass)*10^6,2),nsmall=2)," ppm")
      str7_5 <- ""
    }else{
      tag1 <- h1(strong("3. View results:"), style = "font-size:22.5px;")
    }
    
    if(input$comp == ""){
      url2 <- a("this repository. ", href="https://github.com/briankleiboeker/SLAT", target="_blank")
      url3 <- a("the associated R package.", href = "https://github.com/briankleiboeker/slatR", target="_blank")
      url4 <- a("bkleiboeker@wustl.edu",href = "mailto:bkleiboeker@wustl.edu", target="_blank")
      url5 <- a("Github. ", href="https://github.com/briankleiboeker/SLAT/issues", target="_blank")
      url6 <- a("test dataset",href = "https://github.com/briankleiboeker/SLAT/blob/main/negative_ion_single_peak_test_dataset.csv", target = "_blank")
      tagList(
        p(HTML(paste("","This webapp was designed for use with high-resolution Orbitrap ESI-MS data with the goal of accurately assigining lipid class and structure to user-provided elemental compositions or m/z values, but should also be compatible with most any high-resolution mass spectromerty data provided in a standard format.","", sep = '<br/>'))),
        tags$img(src = "slat_logo.png", height = "300px", width = "500px",alt = "somethingwentwrong",deleteFile=FALSE),
        p(HTML(paste("",paste0("Report problems, suggestions, or other feedback to ",url4, " or open an issue on ",url5),"", sep = '<br/>'))),
        p("Instructions for usage, source code, and a ",url6," are available at ",url2,""," For reproducible, high-throughput assignment of lipid class and structure, check out ",url3),
        p(HTML(paste("","","Input an elemental composition or m/z value to begin",sep = '<br/>')))
      )
    }else{
      url3 <- a("this repository. ", href="https://github.com/briankleiboeker/SLAT", target="_blank")
      url4 <- a("bkleiboeker@wustl.edu",href = "mailto:bkleiboeker@wustl.edu", target="_blank")
      url5 <- a("Github. ", href="https://github.com/briankleiboeker/SLAT/issues", target="_blank")
      url6 <- a("test dataset",href = "https://github.com/briankleiboeker/SLAT/blob/main/negative_ion_single_peak_test_dataset.csv", target = "_blank")
      tagList(
        p(HTML(paste(paste0("Report problems, suggestions, or other feedback to ",url4, " or open an issue on ",url5),"",paste0("Instructions for usage, source code, and a ",url6," are available at ",url3),sep = '<br/>'))),
        tag1,
      p(HTML(paste(" ",str1, str2, str3, str4, str5, str5_5, str6, str7,str7_5,str8,str9,str10,str11,str12,sep = '<br/>')))
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
            HTML(paste0("Import successful, but no column names detected. ", "<b>Assuming that the rightmost non-empty column contains compositions.</b> Filter output by abs(delta) cutoff feature is not available.",ifelse( any(data.full()$structure %in% pls$gen.structure) ,""," Exact structure is unknown for all assigned structures, so exact structure columns are not being shown."),sep = '<br/>'))
          }else{
            HTML(paste0("Import successful, but composition column not detected. ", "<b>Assuming that the rightmost non-empty column contains compositions.</b> ",ifelse("delta" %in% colnames(data.full()),"Delta column detected, so filtering by abs(delta) cutoff feature is enabled. ","Delta column not detected, so filtering by abs(delta) cutoff feature is disabled. "),ifelse( any(data.full()$structure %in% pls$gen.structure) ,""," Exact structure is unknown for all assigned species, so exact structure columns are not being shown."),sep = '<br/>'))
          }
        }else{
          HTML(paste0("Import successful. Detected columns: ", paste( colnames(data.full())[!(colnames(data.full()) %in% c("structure","strain.specific.assignment","species.specific.assignment","Analysis_parameters:")) & !grepl("exact.structure",colnames(data.full()))],collapse = ', '),ifelse("delta" %in% colnames(data.full()),". Delta column detected, so filtering by abs(delta) cutoff feature is enabled. ",". Delta column not detected, so filtering by abs(delta) cutoff feature is disabled. "),ifelse( any(data.full()$structure %in% pls$gen.structure) ,""," Exact structure is unknown for all assigned structures, so exact structure columns are not being shown.")))
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
        df$s <- return_num_elements("S",df$composition)
        df$na <- return_num_elements("Na",df$composition)
        df$cl <- return_num_elements("Cl",df$composition)
        df$li <- return_num_elements("Li",df$composition)
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
        df$s <- return_num_elements("S",df[,ncols])
        df$na <- return_num_elements("Na",df[,ncols])
        df$cl <- return_num_elements("Cl",df[,ncols])
        df$li <- return_num_elements("Li",df[,ncols])
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
        bulkadducts <- c("m.plus.h","m.plus.ammonia","m.plus.sodium","m.minus.h","m.plus.chloride","m.minus.2h","m.plus.lithium")
        bulklois <- input$bulkneutralloi
      }
      
      df$structure <- mapply(assign_species,df$c,df$h,df$o,df$n,df$p,df$s,df$na,df$cl,df$li,input$bulkionmode,rep(list(bulkadducts),nrow(df)),( df$c - ((df$h+df$cl+df$na+df$li)/2) + ((df$n+df$p)/2) +1 ),input$bulkdomain,rep(list(bulklois),nrow(df)),input$bulkmaxdblbnds,input$bulkstrain,SIMPLIFY = T)
      df$structure <- ifelse(is.na(df$structure),"",df$structure)
      df<-df[, !(colnames(df) %in% c("c","h","o","n","p","s","na","cl","li"))]
      
      
      
      
      
      
      if(any(df$structure %in% pls$gen.structure)){
        
        if( (input$bulkstrain != "n/a"  && any(df$structure %in% pls[pls$strain == input$bulkstrain , colnames(pls) == "gen.structure"] )) | (input$bulkspecies != "n/a"  && any(df$structure %in% pls[pls$species == input$bulkspecies , colnames(pls) == "gen.structure"] )) ){
          
          if(input$bulkstrain != "n/a"  && any(df$structure %in% pls[pls$strain == input$bulkstrain , colnames(pls) == "gen.structure"] )){
            pls2 <- pls[pls$strain == input$bulkstrain,]
          }else{ #this bool is true: 'input$bulkspecies != "n/a"  && any(df$structure %in% pls[pls$species == input$bulkspecies , colnames(pls) == "gen.structure"] )'
            pls2 <- pls[pls$species == input$bulkspecies,]
          }

          df$exact.structure.1 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[1],"" )) )
          df$exact.structure.2 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[2],"" )) )
          df$exact.structure.3 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[3],"" )) )
          df$exact.structure.4 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[4],"" )) )
          df$exact.structure.5 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[5],"" )) )
          df$exact.structure.6 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[6],"" )) )
          df$exact.structure.7 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[7],"" )) )
          
          if(input$bulkstrain != "n/a"){
            df$strain.specific.assignment <- ifelse(df$exact.structure.1 != "","Y","N")
          }else{
            df$species.specific.assignment <- ifelse(df$exact.structure.1 != "","Y","N")
          }
          
          
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
          pls2 <- pls[pls$strain == "" & pls$species == "",]
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
      
      if("delta" %in% colnames(df)){
        if(input$returnAll == "returnassignedpeaks"){
          df <- df[ifelse( !is.na(df$delta) & abs(df$delta) <= input$deltacutoff & df$structure != "" ,T,F) ,]
        }else{
          df <- df[ifelse( !is.na(df$delta) & abs(df$delta) <= input$deltacutoff ,T,F) ,]
        }
      }else{
        if(input$returnAll == "returnassignedpeaks"){
          df <- df[ifelse(df$structure != "" ,T,F) ,]
        }
      }
      
      #add analysis details here - recall that this is "assign starting from composition" module
      new_col_number <- (ncol(df) + 1)
      sample_source <- input$bulkdomain 
      strain_species <- ifelse(input$bulkstrain != "n/a" | input$bulkspecies != "n/a",
                               ifelse(input$bulkstrain != "n/a",input$bulkstrain,input$bulkspecies),
                               "n/a")
      ion_mode <- input$bulkionmode
      if(length(input$bulknegionadduct) >= 1 | length(input$bulkposionadduct) >= 1){
        adducts <- ifelse(input$bulkionmode == "neg.ion",
                          paste0(unlist(input$bulknegionadduct),collapse = "; "), #enter here if negion mode
                          ifelse(input$bulkionmode == "pos.ion", #enter here if posionmode or neutralionmode
                                 paste0(unlist(input$bulkposionadduct),collapse = "; "),
                                 "n/a"
                                 )
                          )
      }else{
        adducts <- ifelse(input$bulkionmode == "neutral",
                          "n/a",
                          "none")
      }
      
      if((input$bulkionmode == "neg.ion" & length(input$bulknegionloi) >= 1) | (input$bulkionmode == "pos.ion" & length(input$bulkposionloi) >= 1) | (input$bulkionmode == "neutral" & length(input$bulkneutralloi) >= 1)){
        lipids_of_interest <- ifelse(input$bulkionmode == "neg.ion",
                                     paste0(unlist(input$bulknegionloi),collapse = "; "),
                                     ifelse(input$bulkionmode == "pos.ion",
                                            paste0(unlist(input$bulkposionloi),collapse = "; "),
                                            paste0(unlist(input$bulkneutralloi),collapse = "; ")
                                     )
        )
      }else{
        lipids_of_interest <- "none"
      }
      return_all <- ifelse(input$returnAll == "returnallpeaks","Return all peaks","Only return peaks with a species assigned")
      max_acyl_alkyl_dbl_bonds <- input$bulkmaxdblbnds
      starting_point <- ifelse(input$assignmentmode == "frommz","Start from m/z values","Start from predetermined elemental compositions")
      # below is where it becomes unique to this "start from composition" module
      isotope <- input$bulkisotope
      deltacutoff <- input$deltacutoff
      #add to df
      df[1,new_col_number] <- paste0("Sample source: ",sample_source)
      df[2,new_col_number] <- paste0("Sample strain/species: ",strain_species)
      df[3,new_col_number] <- paste0("Ion mode: ",ion_mode)
      df[4,new_col_number] <- paste0("Ions/adducts searched for: ",adducts)
      df[5,new_col_number] <- paste0("Lipid classes searched for: ",lipids_of_interest)
      df[6,new_col_number] <- return_all
      df[7,new_col_number] <- paste0("Max acyl/alkyl double bonds: ",max_acyl_alkyl_dbl_bonds)
      df[8,new_col_number] <- starting_point
      # below is where it resumes uniqueness to this "start from composition" module
      df[9,new_col_number] <- paste0("Isotope being used: ",isotope)
      df[10,new_col_number] <- paste0("abs(Delta) cutoff: ",deltacutoff)
      # back to shared
      colnames(df)[new_col_number] <- "Analysis_parameters:"
      df[11:nrow(df),new_col_number] <- ""
      
      df
    }else{
      imported$data()
    }
    
  })
  
  output$data <- renderPrint({
    if(!is.null(imported$data()) ){
      if("delta" %in% colnames(data.full())){
        if(input$returnAll == "returnassignedpeaks"){
          data.full()[ifelse( !is.na(data.full()$delta) & abs(data.full()$delta) <= input$deltacutoff & data.full()$structure != "" ,T,F) ,colnames(data.full()) != "Analysis_parameters:"] |> head(n=10)
        }else{
          data.full()[ifelse( !is.na(data.full()$delta) & abs(data.full()$delta) <= input$deltacutoff ,T,F) ,colnames(data.full()) != "Analysis_parameters:"] |> head(n=10)
        }
      }else{
        if(input$returnAll == "returnassignedpeaks"){
          data.full()[ifelse(data.full()$structure != "" ,T,F) ,colnames(data.full()) != "Analysis_parameters:"] |> head(n=10)
        }else{
          data.full()[,colnames(data.full()) != "Analysis_parameters:"] |> head(n=10)
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
      write.csv( data.full(),file,row.names = F)
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
            HTML(paste0("Import successful, but no column names detected. ", "<b>Assuming that the leftmost non-empty column contains m/z values.</b>",ifelse( any(data.full2()$structure %in% pls$gen.structure) ,""," Exact structure is unknown for all assigned structures, so exact structure columns are not being shown."),sep = '<br/>'))
          }else{
            HTML(paste0("Import successful, but m/z column not detected. ", "<b>Assuming that the leftmost non-empty column contains m/z values.</b> ",ifelse( any(data.full2()$structure %in% pls$gen.structure) ,""," Exact structure is unknown for all assigned species, so exact structure columns are not being shown."),sep = '<br/>'))
          }
        }else{
          HTML(paste0("Import successful. Detected columns: ", paste( colnames(data.full2())[ !(colnames(data.full2()) %in% c("structure","composition","theo.mass","rdb","delta","strain.specific.assignment","species.specific.assignment","Analysis_parameters:")) & !grepl("exact.structure",colnames(data.full2()))],collapse = ', '),ifelse( any(data.full2()$structure %in% pls$gen.structure) ,"",". Exact structure is unknown for all assigned structures, so exact structure columns are not being shown.")))
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
            HTML(paste0("Import successful, but no column names detected. ", "<b>Assuming that the leftmost non-empty column contains m/z values.</b>",sep = '<br/>'))
          }else{
            HTML(paste0("Import successful, but m/z column not detected. ", "<b>Assuming that the leftmost non-empty column contains m/z values.</b> ",sep = '<br/>'))
          }
        }else{
          HTML(paste0("Import successful. Detected columns: ", paste( colnames(data.full3())[ !(colnames(data.full3()) %in% c("corr.mz","intensity.plot.col","Recalibration_summary:"))],collapse = ', ')))
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

  observeEvent(event_data("plotly_click", source = "plot1"),{
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
         theme(text=element_text(size=15))
       
       p1 %>% ggplotly(source = "plot1") %>% 
         event_register('plotly_click') %>% config( #displayModeBar = F
           displaylogo = F,
           modeBarButtonsToRemove = c('zoom','pan','select','zoomIn','zoomOut','lasso2d','hoverClosestCartesian', 'hoverCompareCartesian','autoScale','logo'),
           toImageButtonOptions = list(height = NULL,
                                       width = NULL,
                                       format= 'jpeg', # one of png, svg, jpeg, webp
                                       filename= 'recal_plot',
                                       height= 750,
                                       width= 1050,
                                       scale= 1.5)
           )

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
          
          #add recalibration info to df
          newcolno <- ncol(df) + 1
          df[1,newcolno] <- paste0("1-point recalibration performed")
          df[2,newcolno] <- paste0("Observed m/z of selected point: ",as.numeric(alignment_table()$mz[1]),", assigned theoretical mass of selected point: ",as.numeric(alignment_table()$theo.mass[1]))
          df[3,newcolno] <- paste0("The difference of theo.mass and m/z (theo.mass - m/z) was added to the observed m/z of all points in dataset to obtain the final corrected m/z")
          colnames(df)[newcolno] <- "Recalibration_summary:"
          df[4:nrow(df),newcolno] <- ""
        }else{
          yvals <- (as.numeric(alignment_table()$theo.mass) - as.numeric(alignment_table()$mz))
          xvals <- as.numeric(alignment_table()$mz)
          s <- summary(lm(yvals ~ xvals))$coefficients
          yint <- as.numeric(s[rownames(s) == "(Intercept)",colnames(s) == "Estimate"])
          slope <- as.numeric(s[rownames(s) == "xvals",colnames(s) == "Estimate"])
          df$corr.mz <- as.numeric(df$mz) + ((slope*as.numeric(df$mz))+yint)
          
          #add recalibration info to df
          newcolno <- ncol(df) + 1
          df[1,newcolno]  <- paste0(nrow(alignment_table()),"-point recalibration performed")
          #paste0(mapply(function(a,b,c,d) paste0(a,"[",b,"]: ",c," to ",d),this_table2()$element,this_table2()$isotope,this_table2()$minimum,this_table2()$maximum),collapse = "; ")
          df[2,newcolno]  <- paste0("Observed m/z and assigned theoretical mass pairs are as follows in format [m/z, theo.mass]: ",
                                    paste0(mapply(function(a,b) paste0("[",a,", ",b,"]"),
                                                  as.numeric(alignment_table()$mz),as.numeric(alignment_table()$theo.mass)
                                                  ),
                                           collapse = "; ")
                                    )
          df[3,newcolno]  <- paste0("A linear model was fit to calculate error (error = theo.mass - m/z) which was subsequently used to recalibrate all data points")
          df[4,newcolno]  <- paste0("Linear model call was: ",paste0(summary(lm(yvals ~ xvals))$call,collapse = " "),", where \'yvals\' = theo.mass - m/z and \'xvals\' = m/z")
          colnames(df)[newcolno] <- "Recalibration_summary:"
          df[5:nrow(df),newcolno] <- ""
        }
      }else{
        df$corr.mz <- rep(NA,nrow(df))
      }
      #NOTE: cannot remove the "intensity.plot.col" here b/c it is used for plotting -> can only remove it in the download handler
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
    recal_info <- c()
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
          
          #have to redo recal here
          if(nrow(alignment_table()) > 0 && all(!is.na(alignment_table()$theo.mass)) ){
            if(nrow(alignment_table()) == 1){
              df$mz <- as.numeric(df$mz) + (as.numeric(alignment_table()$theo.mass[1]) - as.numeric(alignment_table()$mz[1]))
              
              #add recalibration info to a list to later add to df
              recal_info <- c(rep("",10))
              recal_info[1] <- paste0("1-point recalibration performed")
              recal_info[2] <- paste0("Observed m/z of selected point: ",as.numeric(alignment_table()$mz[1]),", assigned theoretical mass of selected point: ",as.numeric(alignment_table()$theo.mass[1]))
              recal_info[3] <- paste0("The difference of theo.mass and m/z (theo.mass - m/z) was added to the observed m/z of all points in dataset to obtain the final corrected m/z")

            }else{
              yvals <- (as.numeric(alignment_table()$theo.mass) - as.numeric(alignment_table()$mz))
              xvals <- as.numeric(alignment_table()$mz)
              s <- summary(lm(yvals ~ xvals))$coefficients
              yint <- as.numeric(s[rownames(s) == "(Intercept)",colnames(s) == "Estimate"])
              slope <- as.numeric(s[rownames(s) == "xvals",colnames(s) == "Estimate"])
              df$mz <- as.numeric(df$mz) + ((slope*as.numeric(df$mz))+yint)
              
              #add recalibration info to a list to later add to df
              recal_info <- c(rep("",10))
              recal_info[1]  <- paste0(nrow(alignment_table()),"-point recalibration performed")
              recal_info[2]  <- paste0("Observed m/z and assigned theoretical mass pairs are as follows in format [m/z, theo.mass]: ",
                                        paste0(mapply(function(a,b) paste0("[",a,", ",b,"]"),
                                                      as.numeric(alignment_table()$mz),as.numeric(alignment_table()$theo.mass)
                                        ),
                                        collapse = "; ")
              )
              recal_info[3]  <- paste0("A linear model was fit to calculate error (error = theo.mass - m/z) which was subsequently used to recalibrate all data points")
              recal_info[4]  <- paste0("Linear model call was: ",paste0(summary(lm(yvals ~ xvals))$call,collapse = " "),", where \'yvals\' = theo.mass - m/z and \'xvals\' = m/z")
              
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
      
      
      if(input$bulkionmode == "pos.ion"){
        bulkadducts <- input$bulkposionadduct
        bulklois <- input$bulkposionloi
      }
      if(input$bulkionmode == "neg.ion"){
        bulkadducts <- input$bulknegionadduct
        bulklois <- input$bulknegionloi
      }
      if(input$bulkionmode == "neutral"){
        
        # I removed m.minus.2h from adducts list for neutral ion mode on 11/17/23
        # It makes the least sense and would take a lot of work to add...
        
        bulkadducts <- c("m.plus.h","m.plus.ammonia","m.plus.sodium","m.minus.h","m.plus.chloride","m.plus.lithium")
        bulklois <- input$bulkneutralloi
      }
      
      
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
                                           error.ppm = input$bulk.delta.error,
                                           adducts = bulkadducts)
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
                                           error.ppm = input$bulk.delta.error,
                                           adducts = bulkadducts)
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
      df$s <- return_num_elements("S",df$composition) + return_num_elements("S[33]",df$composition) + return_num_elements("S[34]",df$composition) + return_num_elements("S[36]",df$composition)
      df$na <- return_num_elements("Na",df$composition)
      df$cl <- return_num_elements("Cl[35]",df$composition) + return_num_elements("Cl[37]",df$composition)
      df$li <- return_num_elements("Li",df$composition) + return_num_elements("Li[6]",df$composition)
      df$rdb <- ifelse(df$theo.mass != "",df$c - ((df$h + df$cl + df$na + df$li) / 2) + ((df$n+df$p)/2) + 1,"")

      df$structure <- mapply(assign_species,df$c,df$h,df$o,df$n,df$p,df$s,df$na,df$cl,df$li,input$bulkionmode,rep(list(bulkadducts),nrow(df)),( df$c - ((df$h + df$cl + df$na + df$li)/2) + ((df$n+df$p)/2) +1 ),input$bulkdomain,rep(list(bulklois),nrow(df)),input$bulkmaxdblbnds,input$bulkstrain,SIMPLIFY = T)
      df$structure <- ifelse(is.na(df$structure),"",df$structure)
      df<-df[, !(colnames(df) %in% c("c","h","o","n","p","s","na","cl","li"))]

      if(any(df$structure %in% pls$gen.structure)){

        if( (input$bulkstrain != "n/a"  && any(df$structure %in% pls[pls$strain == input$bulkstrain , colnames(pls) == "gen.structure"] )) | (input$bulkspecies != "n/a"  && any(df$structure %in% pls[pls$species == input$bulkspecies , colnames(pls) == "gen.structure"] ))  ){
          
          if(input$bulkstrain != "n/a"  && any(df$structure %in% pls[pls$strain == input$bulkstrain , colnames(pls) == "gen.structure"] )){
            pls2 <- pls[pls$strain == input$bulkstrain,]
          }else{
            pls2 <- pls[pls$species == input$bulkspecies,]
          }
          
          
          df$exact.structure.1 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[1],"" )) )
          df$exact.structure.2 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[2],"" )) )
          df$exact.structure.3 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[3],"" )) )
          df$exact.structure.4 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[4],"" )) )
          df$exact.structure.5 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[5],"" )) )
          df$exact.structure.6 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[6],"" )) )
          df$exact.structure.7 <- unlist(lapply(df$structure,function(x) ifelse(x %in% pls2$gen.structure,unname(unlist(pls2[match(x,pls2$gen.structure),2:ncol(pls2)]))[7],"" )) )
          
          if(input$bulkstrain != "n/a"){
            df$strain.specific.assignment <- ifelse(df$exact.structure.1 != "","Y","N")
          }else{
            df$species.specific.assignment <- ifelse(df$exact.structure.1 != "","Y","N")
          }
          
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
          pls2 <- pls[pls$strain == "" & pls$species == "",]
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
      
      if(input$returnAll == "returnassignedpeaks"){
        df <- df[df$structure != "",]
      }

      #add analysis details here - recall that this is "assign starting from m/z" module
      
      new_col_number <- (ncol(df) + 1)
      sample_source <- input$bulkdomain
      strain_species <- ifelse(input$bulkstrain != "n/a" | input$bulkspecies != "n/a",
                               ifelse(input$bulkstrain != "n/a",input$bulkstrain,input$bulkspecies),
                               "n/a")
      ion_mode <- input$bulkionmode
      if(length(input$bulknegionadduct) >= 1 | length(input$bulkposionadduct) >= 1){
        #enter here if there is a posionadduct or a negionadduct selected
        adducts <- ifelse(input$bulkionmode == "neg.ion",
                          paste0(unlist(input$bulknegionadduct),collapse = "; "), #enter here if negion mode
                          ifelse(input$bulkionmode == "pos.ion", #enter here if posionmode or neutralionmode
                                 paste0(unlist(input$bulkposionadduct),collapse = "; "),
                                 "n/a"
                          )
        )
      }else{
        adducts <- ifelse(input$bulkionmode == "neutral",
                          "n/a",
                          "none")
      }
      
      if((input$bulkionmode == "neg.ion" & length(input$bulknegionloi) >= 1) | (input$bulkionmode == "pos.ion" & length(input$bulkposionloi) >= 1) | (input$bulkionmode == "neutral" & length(input$bulkneutralloi) >= 1)){
        lipids_of_interest <- ifelse(input$bulkionmode == "neg.ion",
                                     paste0(unlist(input$bulknegionloi),collapse = "; "),
                                     ifelse(input$bulkionmode == "pos.ion",
                                            paste0(unlist(input$bulkposionloi),collapse = "; "),
                                            paste0(unlist(input$bulkneutralloi),collapse = "; ")
                                     )
        )
      }else{
        lipids_of_interest <- "none"
      }
      return_all <- ifelse(input$returnAll == "returnallpeaks","Return all peaks","Only return peaks with a species assigned")
      max_acyl_alkyl_dbl_bonds <- input$bulkmaxdblbnds
      starting_point <- ifelse(input$assignmentmode == "frommz","Start from m/z values","Start from predetermined elemental compositions")
      # below is where it becomes unique to this "start from m/z" module
      rdbrule_var <- input$bulkrdbrule
      rdbrange_var <-  paste0(as.numeric(input$bulkrdbrange[1])," to ",as.numeric(input$bulkrdbrange[2]))
      delta_error_ppm <- input$bulk.delta.error
      element_table <- paste0(mapply(function(a,b,c,d) paste0(a,"[",b,"]: ",c," to ",d),this_table2()$element,this_table2()$isotope,this_table2()$minimum,this_table2()$maximum),collapse = "; ")
      #add to df
      df[1,new_col_number] <- paste0("Sample source: ",sample_source)
      df[2,new_col_number] <- paste0("Sample strain/species: ",strain_species)
      df[3,new_col_number] <- paste0("Ion mode: ",ion_mode)
      df[4,new_col_number] <- paste0("Ions/adducts searched for: ",adducts)
      df[5,new_col_number] <- paste0("Lipid classes searched for: ",lipids_of_interest)
      df[6,new_col_number] <- return_all
      df[7,new_col_number] <- paste0("Max acyl/alkyl double bonds: ",max_acyl_alkyl_dbl_bonds)
      df[8,new_col_number] <- starting_point
      # below is where it resumes uniqueness to this "start from composition" module
      df[9,new_col_number] <- paste0("RDB restriction rule: ",rdbrule_var)
      df[10,new_col_number] <- paste0("RDB restriction range: ",rdbrange_var)
      df[11,new_col_number] <- paste0("Mass tolerance (ppm): ",delta_error_ppm)
      df[12,new_col_number] <- paste0("Element table for determining composition from m/z: ",element_table)
      # back to shared
      colnames(df)[new_col_number] <- "Analysis_parameters:"
      df[13:nrow(df),new_col_number] <- ""
      
      if(length(recal_info)>0){
        df[1:length(recal_info),new_col_number+1] <- recal_info
        colnames(df)[new_col_number+1] <- "Recalibration_summary:"
        df[(length(recal_info)+1):nrow(df),new_col_number+1] <- ""
      }
      
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
        data.full2()[data.full2()$structure != "",!(colnames(data.full2()) %in% c("Analysis_parameters:","Recalibration_summary:"))] |> head(n = 10)
        }else{
          data.full2()[,!(colnames(data.full2()) %in% c("Analysis_parameters:","Recalibration_summary:"))] |> head(n = 10)
          }
    }
      
  })
  
  output$clip <- renderUI({
    if(!is.null(data.full3())){
      rclipButton(inputId = "clipbtn",
                  label = "Copy aligned data to clipboard",
                  clipText = readr::format_tsv(suppressWarnings(dplyr::rename(dplyr::select(data.full3()[,!(colnames(data.full3()) %in% c("intensity.plot.col"))],!c("mz")),"mz" = "corr.mz") )),
                  icon = icon("clipboard"),
                  style="color: #FFFFFF; background-color: #7DA8F2; border-color: #7DA8F2")

    }else{
      rclipButton(inputId = "clipbtn",
                  label = "Copy aligned data to clipboard", 
                  clipText = "something went wrong",
                  icon = icon("clipboard"),
                  style="color: #FFFFFF; background-color: #7DA8F2; border-color: #7DA8F2")
      
    }

  })
  
  output$download2 <- downloadHandler(
    filename = function(){
      paste0(imported2$name(),".csv")
    },
    content = function(file){
      write.csv(data.full2(),file,row.names = F)
    }
  )
  
  #testdataset handlers:
  output$testdataset1 <- downloadHandler(
    filename = function(){
      "eukaryote_negative_ion_testdata_2.csv"
    },
    content = function(file){
      file.copy("negative_ion_testdata_2.csv", file)
    }
  )
  output$testdataset2 <- downloadHandler(
    filename = function(){
      "eukaryote_negative_ion_testdata_1.csv"
    },
    content = function(file){
      file.copy("negative_ion_testdata_1.csv", file)
    }
  )
  
  output$testdataset3 <- downloadHandler(
    filename = function(){
      "eukaryote_negative_ion_testdata_1.csv"
    },
    content = function(file){
      file.copy("negative_ion_testdata_1.csv", file)
    }
  )
  
  #recalibration tab data download handler
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
