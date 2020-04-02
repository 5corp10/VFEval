# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Initialize package (Run first)
#'
#' Import dependent packages, import patient data and point mapping, and set constants
#' @export
pkgInit = function(data_format = 2)
{
  #import libraries
  library("tidyverse")
  library("dplyr")
  library("ggplot2")
  library("readxl")
  library("devtools")
  library("roxygen2")
  library("cowplot")
  library("UpSetR")
  library("venneuler")
  library("readxl")
  library("reshape2")
  library("scales")
  library("spatstat.utils")
  library("binom")
  library("gridGraphics")

  if(data_format == 1){
    # specify excel file
    excel.file = file.path("Auk_VFs.xlsx")

    # import all data
    assign("df.vf_data", read_xlsx(excel.file, sheet = 1), envir = .GlobalEnv)

    # import point mapping
    assign("df.pt_mapping", read_xlsx(excel.file, sheet = 2), envir = .GlobalEnv)

    # set indexing constants for excising VF data from df.vf_data
    assign("VF_V_OFST", 2, envir = .GlobalEnv)
    assign("VF_V_SIZE", 1, envir = .GlobalEnv)
    assign("T_H_OFST", 56, envir = .GlobalEnv)
    assign("VF_H_SIZE", 76, envir = .GlobalEnv)
    assign("TD_H_OFST", T_H_OFST + VF_H_SIZE + 1, envir = .GlobalEnv)
    assign("TDP_H_OFST", TD_H_OFST + VF_H_SIZE + 1, envir = .GlobalEnv)
    assign("PD_H_OFST", TDP_H_OFST + VF_H_SIZE + 1, envir = .GlobalEnv)
    assign("PDP_H_OFST", PD_H_OFST + VF_H_SIZE + 1, envir = .GlobalEnv)

    assign("FIXLOSS_H_OFST", 33, envir = .GlobalEnv)
    assign("FPOSERR_H_OFST", 34, envir = .GlobalEnv)
    assign("FNEGERR_H_OFST", 36, envir = .GlobalEnv)

    assign("GHT_H_OFST", 51, envir = .GlobalEnv)
    assign("VFI_H_OFST", 50, envir = .GlobalEnv)
    assign("MD_H_OFST", 47, envir = .GlobalEnv)
    assign("PSD_H_OFST", 49, envir = .GlobalEnv)

    # read patient number
    assign("NUM_PAT", nrow(df.vf_data) - VF_V_OFST, envir = .GlobalEnv)
  }
  else if(data_format == 2){
    # set indexing constants for excising VF data from df.vf_data
    assign("VF_V_OFST", 0, envir = .GlobalEnv)
    assign("VF_V_SIZE", 1, envir = .GlobalEnv)
    assign("T_H_OFST", 50, envir = .GlobalEnv)
    assign("VF_H_SIZE", 76, envir = .GlobalEnv)
    assign("TD_H_OFST", T_H_OFST + VF_H_SIZE + 1, envir = .GlobalEnv)
    assign("TDP_H_OFST", TD_H_OFST + VF_H_SIZE + 1, envir = .GlobalEnv)
    assign("PD_H_OFST", TDP_H_OFST + VF_H_SIZE + 1, envir = .GlobalEnv)
    assign("PDP_H_OFST", PD_H_OFST + VF_H_SIZE + 1, envir = .GlobalEnv)

    assign("FIXLOSS_H_OFST", 27, envir = .GlobalEnv)
    assign("FPOSERR_H_OFST", 30, envir = .GlobalEnv)
    assign("FNEGERR_H_OFST", 34, envir = .GlobalEnv)

    assign("GHT_H_OFST", 45, envir = .GlobalEnv)
    assign("VFI_H_OFST", 44, envir = .GlobalEnv)
    assign("MD_H_OFST", 41, envir = .GlobalEnv)
    assign("PSD_H_OFST", 43, envir = .GlobalEnv)

    # specify excel file and import data
    print("Importing and parsing VF data...")
    excel.file = file.path("../Initial_March_trimmed/VF2016_ID_24-2.xlsx")
    df.Vf_2016 = read_xlsx(excel.file, sheet = 1, range = cell_rows(c(5, NA)), col_names = F)
    excel.file = file.path("../Initial_March_trimmed/VF2017_ID_24-2.xlsx")
    df.Vf_2017 = read_xlsx(excel.file, sheet = 1, range = cell_rows(c(5, NA)), col_names = F)
    excel.file = file.path("../Initial_March_trimmed/VF2018_ID_24-2.xlsx")
    df.Vf_2018 = read_xlsx(excel.file, sheet = 1, range = cell_rows(c(5, NA)),col_names = F)

    # make one large data frame of all Vf data
    assign("df.vf_data", rbind(df.Vf_2016, df.Vf_2017, df.Vf_2018), envir = .GlobalEnv)

    # apply exclusion criteria
    print("Applying exclusion criteria...")
    excel.file = file.path("../ExclusionsClean.csv")
    df.exclusions = read.csv(excel.file)
    list.excl = vector()
    for(excl_i in 1:nrow(df.exclusions)){
      #print(excl_i)
      for(vf_i in which(df.vf_data[,1] == df.exclusions[excl_i,"StudyID"])){
        #print(vf_i)
        eye_excl = df.exclusions[excl_i,"ExcludeEyes"]
        eye_vf = df.vf_data[vf_i,4]
        if(((eye_excl == "OD") && (eye_vf == "R")) || ((eye_excl == "OS") && (eye_vf == "L")) || (eye_excl == "OU")){
          list.excl = c(list.excl,vf_i)
        }
      }
    }
    df.vf_data = df.vf_data[-list.excl, ]

    # read patient number
    assign("NUM_PAT", nrow(df.vf_data) - VF_V_OFST, envir = .GlobalEnv)

    #add extra columns at the end because they are not imported since they are blank
    assign("df.vf_data", cbind(df.vf_data, data.frame(matrix(NA, nrow = NUM_PAT, ncol = 5))), envir = .GlobalEnv)


    # import rnfl, mrw and gcl data
    print("Importing rnfl, mrw and gcl data...")
    excel.file = file.path("../Initial_March_trimmed/RNFL_Q15_mm1_D3436.csv")
    assign("df.rnfl_data", read.csv(excel.file), envir = .GlobalEnv)
    excel.file = file.path("../Initial_March_trimmed/MRW_Q15_Disp.csv")
    assign("df.mrw_data", read.csv(excel.file), envir = .GlobalEnv)
    excel.file = file.path("../Initial_March_trimmed/Macula_Q15_GCL.csv")
    assign("df.gcl_data", read.csv(excel.file), envir = .GlobalEnv)

    #format GCL data - assumed it is ordered
    df.gcl = data.frame(matrix(NA, nrow = 0, ncol = 11))
    for(pat in seq(1,nrow(df.gcl_data),7)){
      df.gcl = rbind(df.gcl, c(df.gcl_data[pat,"StudyID"], df.gcl_data[pat,"Eye"],
                               as.integer(as.POSIXct(df.gcl_data[pat,"ExamDate"])), df.gcl_data[pat,"Quality"],
                               df.gcl_data[pat,"abnormalSector"], df.gcl_data[(pat + 4),"abnormalSector"],
                               df.gcl_data[(pat + 6),"abnormalSector"], df.gcl_data[(pat + 3),"abnormalSector"],
                               df.gcl_data[(pat + 1),"abnormalSector"], df.gcl_data[(pat + 5),"abnormalSector"],
                               df.gcl_data[(pat + 2),"abnormalSector"]))

    }

    colnames(df.gcl) = c("StudyID", "Eye", "ExamDate", "Quality",
                         "GCL.ABN.Global", "GCL.ABN.Sup", "GCL.ABN.TS", "GCL.ABN.NS",
                         "GCL.ABN.Inf", "GCL.ABN.NI", "GCL.ABN.TI")
    assign("df.gcl_data", df.gcl, envir = .GlobalEnv)

    # import helper tables
    print("Importing helper tables...")
    excel.file = file.path("helper_tables.xlsx")
    assign("df.coord", data.frame(read_xlsx(excel.file, sheet = 1, col_names = T)), envir = .GlobalEnv)
    assign("df.pt_mapping", data.frame(read_xlsx(excel.file, sheet = 2, col_names = T)), envir = .GlobalEnv)

    print("Init complete.")
  }
}

#' Partial STATPAC SFA
#'
#' Print reliability and global indices, criteria results, and VF maps (TD, TDP, PD, PDP) for a given patient
#' @param pat_id Patient position in the list on the excel spreadsheet
#' @return grid_plot The grid plot of the STATPAC SFA analysis.
#' @export
printVfMap = function(pat_id)
{
  # extract reliability indices
  fix_loss1 = df.vf_data[VF_V_OFST+pat_id,FIXLOSS_H_OFST]
  fix_loss2 = df.vf_data[VF_V_OFST+pat_id,FIXLOSS_H_OFST+1]
  f_pos_err = percent(df.vf_data[VF_V_OFST+pat_id,FPOSERR_H_OFST])
  f_neg_err = percent(df.vf_data[VF_V_OFST+pat_id,FNEGERR_H_OFST])

  list1 = c(paste0("Fixation losses: ", fix_loss1, "/", fix_loss2), paste0("False POS errors: ", f_pos_err), paste0("False NEG errors: ", f_neg_err))

  # extract global indices
  ght = df.vf_data[VF_V_OFST+pat_id,GHT_H_OFST]
  vfi = percent(df.vf_data[VF_V_OFST+pat_id,VFI_H_OFST])
  md = df.vf_data[VF_V_OFST+pat_id,MD_H_OFST-1]
  if(!is.na(df.vf_data[VF_V_OFST+pat_id,MD_H_OFST]))
    md_p = percent(df.vf_data[VF_V_OFST+pat_id,MD_H_OFST])
  else
    md_p = "NA"
  psd = df.vf_data[VF_V_OFST+pat_id,(PSD_H_OFST-1)]
  if(!is.na(df.vf_data[VF_V_OFST+pat_id,PSD_H_OFST]))
    psd_p = percent(df.vf_data[VF_V_OFST+pat_id,PSD_H_OFST])
  else
    psd_p = "NA"

  list2 = c(paste0("GHT: ", ght), paste0("VFI: ", vfi), paste0("MD: ", md, " dB (<", md_p, ")"), paste0("PSD: ", psd, " dB (<", psd_p, ")"))

  indices = ggdraw() + draw_text(list1, x = 0.2, y = c(0.9, 0.8, 0.7), size = 12, hjust = 0) + draw_text(list2, x = 0.2, y = c(0.55, 0.45, 0.35, 0.25), size = 12, hjust = 0)

  # apply criteria and generate result
  list3 = c(paste0("GHT: ", checkGhtCriteria(pat_id)), paste0("FOST: ", checkFostCriteria(pat_id)),
            paste0("HAP2: ", checkMhpaCriteria(pat_id)), paste0("LOGTS: ", checkLogtsCriteria(pat_id)),
            paste0("UKGTS: ", checkUkgtsCriteria(pat_id)))
  results = ggdraw() + draw_text(list3, x = 0.2, y = c(0.8, 0.7, 0.6, 0.5, 0.4), size = 12, hjust = 0)

  # generate plot of threshold
  df.T = data.frame(df.coord, t(df.vf_data[VF_V_OFST+pat_id, T_H_OFST:(T_H_OFST+VF_H_SIZE-1)]))
  colnames(df.T) = c("x", "y", "T")
  plot.T = ggplot(df.T, aes(x, y, label = T)) + geom_point(alpha=0) + geom_text(aes(label=T),hjust=0.5, vjust=0.5) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

  # generate plot of total deviation
  df.TD = data.frame(df.coord, t(df.vf_data[VF_V_OFST+pat_id, TD_H_OFST:(TD_H_OFST+VF_H_SIZE-1)]))
  colnames(df.TD) = c("x", "y", "TD")
  plot.TD = ggplot(df.TD, aes(x, y, label = TD)) + geom_point(alpha=0) + geom_text(aes(label=TD),hjust=0.5, vjust=0.5) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

  # generate plot of pattern deviation
  df.PD = data.frame(df.coord, t(df.vf_data[VF_V_OFST+pat_id, PD_H_OFST:(PD_H_OFST+VF_H_SIZE-1)]))
  colnames(df.PD) = c("x", "y", "PD")
  plot.PD = ggplot(df.PD, aes(x, y, label = PD)) + geom_point(alpha=0) + geom_text(aes(label=PD),hjust=0.5, vjust=0.5) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

  # generate plot of total deviation probability
  df.TDp = data.frame(df.coord, t(df.vf_data[VF_V_OFST+pat_id, TDP_H_OFST:(TDP_H_OFST+VF_H_SIZE-1)]))
  colnames(df.TDp) = c("x", "y", "TDp")
  plot.TDp = ggplot(df.TDp, aes(x, y, label = TDp)) + geom_point(alpha=0) + geom_text(aes(label=TDp),hjust=0.5, vjust=0.5) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

  # generate plot of pattern deviation probability
  df.PDp = data.frame(df.coord, t(df.vf_data[VF_V_OFST+pat_id, PDP_H_OFST:(PDP_H_OFST+VF_H_SIZE-1)]))
  colnames(df.PDp) = c("x", "y", "PDp")
  plot.PDp = ggplot(df.PDp, aes(x, y, label = PDp)) + geom_point(alpha=0) + geom_text(aes(label=PDp),hjust=0.5, vjust=0.5) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)

  return(plot_grid(indices, results, plot.TD, plot.PD, plot.TDp, plot.PDp, ncol = 2, labels = c("", "", "TD", "PD", "TDp", "PDp")))
}

#' Euler diagram
#'
#' Print euler diagram (i.e. weighted venn diagram) for the overlapping criteria results groups
#' @export
printEulerDiag = function(var_category = "all", df.results = df.best_match)
{
  o = 0
  f = m = u = g = l = 0
  f_m = f_u = f_g = f_l = m_u = m_g = m_l = u_g = u_l = g_l = 0
  f_m_u = f_m_g = f_m_l = f_u_g = f_u_l = f_g_l = m_u_g = m_u_l = m_g_l = u_g_l = 0
  f_m_u_g = f_m_u_l = f_m_g_l = f_u_g_l = m_u_g_l = 0
  f_m_u_g_l = 0

  for(pat in 1:nrow(df.results)){
    #print(pat)
    fost = df.results[pat,"FOST"]
    HAP2 = df.results[pat,"HAP2"]
    ukgts = df.results[pat,"UKGTS"]
    ght = df.results[pat,"GHT"]
    logts = df.results[pat,"LOGTS"]

    if((fost == T) && (HAP2 == T) && (ukgts == T) && (ght == T) && (logts == T))
      f_m_u_g_l = f_m_u_g_l + 1
    else if((fost == T) && (HAP2 == T) && (ukgts == T) && (ght == T))
      f_m_u_g = f_m_u_g + 1
    else if((fost == T) && (HAP2 == T) && (ukgts == T) && (logts == T))
      f_m_u_l = f_m_u_l + 1
    else if((fost == T) && (HAP2 == T) && (ght == T) && (logts == T))
      f_m_g_l = f_m_g_l + 1
    else if((fost == T) && (ukgts == T) && (ght == T) && (logts == T))
      f_u_g_l = f_u_g_l + 1
    else if((HAP2 == T) && (ukgts == T) && (ght == T) && (logts == T))
      m_u_g_l = m_u_g_l + 1
    else if ((fost == T) && (HAP2 == T) && (ukgts == T))
      f_m_u = f_m_u + 1
    else if ((fost == T) && (HAP2 == T) && (ght == T))
      f_m_g = f_m_g + 1
    else if ((fost == T) && (HAP2 == T) && (logts == T))
      f_m_l = f_m_l + 1
    else if ((fost == T) && (ukgts == T) && (ght == T))
      f_u_g = f_u_g + 1
    else if ((fost == T) && (ukgts == T) && (logts == T))
      f_u_l = f_u_l + 1
    else if ((fost == T) && (ght == T) && (logts == T))
      f_g_l = f_g_l + 1
    else if ((HAP2 == T) && (ukgts == T) && (ght == T))
      m_u_g = m_u_g + 1
    else if ((HAP2 == T) && (ukgts == T) && (logts == T))
      m_u_l = m_u_l + 1
    else if ((HAP2 == T) && (ght == T) && (logts == T))
      m_g_l = m_g_l + 1
    else if ((ukgts == T) && (ght == T) && (logts == T))
      u_g_l = u_g_l + 1
    else if ((fost == T) && (HAP2 == T))
      f_m = f_m + 1
    else if ((fost == T) && (ukgts == T))
      f_u = f_u + 1
    else if ((fost == T) && (ght == T))
      f_g = f_g + 1
    else if ((fost == T) && (logts == T))
      f_l = f_l + 1
    else if ((HAP2 == T) && (ukgts == T))
      m_u = m_u + 1
    else if ((HAP2 == T) && (ght == T))
      m_g = m_g + 1
    else if ((HAP2 == T) && (logts == T))
      m_l = m_l + 1
    else if ((ukgts == T) && (ght == T))
      u_g = u_g + 1
    else if ((ukgts == T) && (logts == T))
      u_l = u_l + 1
    else if ((ght == T) && (logts == T))
      g_l = g_l + 1
    else if ((fost == T))
      f = f + 1
    else if ((HAP2 == T))
      m = m + 1
    else if ((ukgts == T))
      u = u + 1
    else if ((ght == T))
      g = g + 1
    else if ((logts == T))
      l = l + 1
    else
      o = o + 1
  }
  areas = c("FOST"=f, "HAP2"=m, "UKGTS"=u, "GHT"=g, "LOGTS"=l,
            "FOST&HAP2"=f_m, "FOST&UKGTS"=f_u, "FOST&GHT"=f_g, "FOST&LOGTS"=f_l, "HAP2&UKGTS"=m_u,
            "HAP2&GHT"=m_g, "HAP2&LOGTS"=m_l, "UKGTS&GHT"=u_g, "UKGTS&LOGTS"=u_l, "GHT&LOGTS"=g_l,
            "FOST&HAP2&UKGTS"=f_m_u, "FOST&HAP2&GHT"=f_m_g, "FOST&HAP2&LOGTS"=f_m_l, "FOST&UKGTS&GHT"=f_u_g, "FOST&UKGTS&LOGTS"=f_u_l,
            "FOST&GHT&LOGTS"=f_g_l, "HAP2&UKGTS&GHT"=m_u_g, "HAP2&GHT&LOGTS"=m_g_l, "HAP2&UKGTS&LOGTS"=m_u_l, "UKGTS&GHT&LOGTS"=u_g_l,
            "FOST&HAP2&UKGTS&GHT"=f_m_u_g, "FOST&HAP2&UKGTS&LOGTS"=f_m_u_l, "FOST&HAP2&GHT&LOGTS"=f_m_g_l, "FOST&UKGTS&GHT&LOGTS"=f_u_g_l, "HAP2&UKGTS&GHT&LOGTS"=m_u_g_l,
            "FOST&HAP2&UKGTS&GHT&LOGTS"=f_m_u_g_l)

  #areas = c("F1"=0,"M1"=0,"U1"=0,"G1"=0,"F1&M1"=0,"F1&U1"=0,"F1&G1"=0,"M1&U1"=0,"M1&G1"=0,"U1&G1"=0,"F1&M1&U1"=0,"F1&M1&G1"=0,"F1&U1&G1"=0,"M1&U1&G1"=0,"F1&M1&U1&G1"=0,
  #          "F2"=0,"M2"=0,"U2"=0,"G2"=0,"F2&M2"=0,"F2&U2"=0,"F2&G2"=0,"M2&U2"=0,"M2&G2"=0,"U2&G2"=0,"F2&M2&U2"=0,"F2&M2&G2"=0,"F2&U2&G2"=0,"M2&U2&G2"=0,"F2&M2&U2&G2"=0,
  #          "F3"=0,"M3"=0,"U3"=0,"G3"=0,"F3&M3"=0,"F3&U3"=0,"F3&G3"=0,"M3&U3"=0,"M3&G3"=0,"U3&G3"=0,"F3&M3&U3"=0,"F3&M3&G3"=0,"F3&U3&G3"=0,"M3&U3&G3"=0,"F3&M3&U3&G3"=0,
  #          "F4"=0,"M4"=0,"U4"=0,"G4"=0,"F4&M4"=0,"F4&U4"=0,"F4&G4"=0,"M4&U4"=0,"M4&G4"=0,"U4&G4"=0,"F4&M4&U4"=0,"F4&M4&G4"=0,"F4&U4&G4"=0,"M4&U4&G4"=0,"F4&M4&U4&G4"=0,
  #          "F5"=0,"M5"=0,"U5"=0,"G5"=0,"F5&M5"=0,"F5&U5"=0,"F5&G5"=0,"M5&U5"=0,"M5&G5"=0,"U5&G5"=0,"F5&M5&U5"=0,"F5&M5&G5"=0,"F5&U5&G5"=0,"M5&U5&G5"=0,"F5&M5&U5&G5"=0,
  #          "F6"=0,"M6"=0,"U6"=0,"G6"=0,"F6&M6"=0,"F6&U6"=0,"F6&G6"=0,"M6&U6"=0,"M6&G6"=0,"U6&G6"=0,"F6&M6&U6"=0,"F6&M6&G6"=0,"F6&U6&G6"=0,"M6&U6&G6"=0,"F6&M6&U6&G6"=0)

  df.interx[,var_category] <<- c(f, m, u, g, l,
                                 f_m, f_u, f_g, f_l, m_u, m_g, m_l, u_g, u_l, g_l,
                                 f_m_u, f_m_g, f_m_l, f_u_g, f_u_l, f_g_l, m_u_g, m_g_l, m_u_l, u_g_l,
                                 f_m_u_g, f_m_u_l, f_m_g_l, f_u_g_l, m_u_g_l,
                                 f_m_u_g_l,
                                 o)

  v = venneuler(areas)
  v$labels <- c("", "", "", "")
  plot(v)
  text(.5, 0.1, var_category, cex=2)

  #text(.4, 0.9, paste0("F&G&M&U(grey) = ", f_m_u_g), cex=1)
  #text(.4, 0.85, paste0("M&U(green) = ", m_u), cex=1)
  #text(.6, 0.9, paste0("M(yellow) = ", m), cex=1)
  #text(.6, 0.85, paste0("U(blue) = ", u), cex=1)
}

#' Intersetion bar graph
#'
#' Print intersection bar graph to compare overlapping criteria results groups
#' @export
printIntersectionBarGraph = function(df.results = df.best_match)
{
  fost = double(1)
  HAP2 = double(1)
  ukgts = double(1)
  ght = double(1)

  i = 1
  j = 1
  k = 1
  l = 1

  for(pat in 1:NUM_PAT){
    if(df.results[pat,"FOST"] == T){
      fost[i] = pat
      i = i + 1
    }
    if(df.results[pat,"HAP2"] == T){
      HAP2[j] = pat
      j = j + 1
    }
    if(df.results[pat,"UKGTS"] == T){
      ukgts[k] = pat
      k = k + 1
    }
    if(df.results[pat,"GHT"] == T){
      ght[l] = pat
      l = l + 1
    }
  }

  list_input = list("FOST" = fost, "HAP2" = HAP2, "UKGTS" = ukgts, "GHT" = ght)
  upset(fromList(list_input), sets = c("GHT", "UKGTS", "HAP2", "FOST"), keep.order = T, empty.intersections = T, order.by = "freq", text.scale = 2)
}

#' Clustered histogram for MD and criteria results
#'
#' Print histogram where criteria (Foster, HAP2, UKGTS, and GHT) results are clustered within their MD percentiles (10, 5, 2, 1, and 0.5 %)
#' @export
printClusteredHist = function(df.results = df.best_match, x_var = "os")
{
  df.results = df.best_match
  names_col = c("HAP2", "UKGTS", "GHT","FOST", "LOGTS", "EAGLE", "AGIS")

  if(x_var == "os"){
    # create data frame to house graph variables
    df.results_os_cluster = data.frame(matrix(0, nrow=4, ncol=length(names_col)+2))
    colnames(df.results_os_cluster) = c("OCT.Score", names_col, "Total")
    df.results_os_cluster[,1] = c("0", "1-3", "4-5","6")
    md.db_0 = md.db_13 = md.db_45 = md.db_6 = double()

    # record counts of OCT scores in the cells
    for(pat in 1:nrow(df.results)){
      os = df.results[pat,"OCT.Score"]
      #print(pat)

      if(os == 0){
        row = 1
        md.db_0 = c(md.db_0, round(as.double(as.character(df.results[pat,"MD.db"])), 2))
      }
      else if(inside.range(os, c(1,3))){
        row = 2
        md.db_13 = c(md.db_13, round(as.double(as.character(df.results[pat,"MD.db"])), 2))
      }
      else if(inside.range(os, c(4,5))){
        row = 3
        md.db_45 = c(md.db_45, round(as.double(as.character(df.results[pat,"MD.db"])), 2))
      }
      else if(os == 6){
        row = 4
        md.db_6 = c(md.db_6, round(as.double(as.character(df.results[pat,"MD.db"])), 2))
      }

      # increment count of OCT bin
      df.results_os_cluster[row,"Total"] = df.results_os_cluster[row,"Total"] + 1

      # increment count of OCT-VFcriterion cell
      for(col in names_col){
        #print(pat)
        if(df.results[pat,col] == T)
          df.results_os_cluster[row,col] = df.results_os_cluster[row,col] + 1

      }
    }
    assign("df.results_os_cluster", df.results_os_cluster, envir = .GlobalEnv)

    # confidence intervals
    CI1 = binom.confint(x=df.results_os_cluster[,2:5], n=df.results_os_cluster[,"Total"], methods="wilson")
    CI2 = binom.confint(x=df.results_os_cluster[,6], n=df.results_os_cluster[,"Total"], methods="wilson")

    # modify name of each bin
    md.db = list(md.db_0, md.db_13, md.db_45, md.db_6)
    for(row in 1:nrow(df.results_os_cluster)){
      df.results_os_cluster[row,1] = paste0(df.results_os_cluster[row,1], "\nmed=", round(median(md.db[[row]]),1), "dB", "\nN=", df.results_os_cluster[row,"Total"])
    }

    # generate data frame to be graphed
    df.results_graph = cbind(OCT.Score=df.results_os_cluster[,1], CI1[7:10], CI2[,4])
    colnames(df.results_graph)[2:6] = names_col[1:5]
    df.results_graph = within(df.results_graph,  OCT.Score <- factor(OCT.Score, levels=OCT.Score))
    print(df.results_graph)

    # melt
    df.melted = melt(df.results_graph, variable.name = "criterion", value.name = "Hit.Rate")
    # round hit rates to 2 sig figs
    df.melted[,"Hit.Rate"] = round(df.melted[,"Hit.Rate"], digits=2)
    # add CI columns
    df.melted = cbind(df.melted,
                      lower.CI=c(CI1[,"lower.HAP2"], CI1[,"lower.UKGTS"], CI1[,"lower.GHT"], CI1[,"lower.FOST"], CI2[,"lower"]),
                      upper.CI=c(CI1[,"upper.HAP2"], CI1[,"upper.UKGTS"], CI1[,"upper.GHT"], CI1[,"upper.FOST"], CI2[,"upper"]))
    print(df.melted)
    plot.hist = ggplot(df.melted, aes(x=OCT.Score, y=Hit.Rate, fill=criterion)) +
      geom_linerange(position=position_dodge(0.5), aes(ymin=lower.CI, ymax=upper.CI), size=0.8) +
      geom_point(position = position_dodge(0.5), stat = "identity", aes(fill = criterion), size = 5, shape = 21, colour = "black", size = 5, stroke = 1.3) +
      scale_fill_manual(values = c(GHT = "#bbbcbe", FOST = "#ffffb1", HAP2 = "#ffb1b1", UKGTS = "#b1e6fa", LOGTS = "white")) +
      #geom_text(aes(label = Repeat.Error.Rate, group = criterion), size=6, hjust=0.5, vjust=3, position=position_dodge(0.9)) +
      theme_bw(base_size = 22) +
      ylim(0,1)
      #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

  }
  else if(x_var == "md"){
    # create data frame to house graph variables
    df.results_md_cluster = data.frame(matrix(0, nrow=4, ncol=length(names_col)+2))
    colnames(df.results_md_cluster) = c("MD", names_col, "Total")
    df.results_md_cluster[,1] = c(">10%", "2-10%", "0.5-2%", "<0.5%")
    md.db_10 = md.db_5 = md.db_1 = md.db_05 = double()

    for(pat in 1:nrow(df.results)){
      md = df.results[pat,"MD.pval"]
      #print(md)

      if(is.na(md)){
        row = 1
        md.db_10 = c(md.db_10, round(as.double(as.character(df.results[pat,"MD.db"])), 2))
      }
      else if(md == 0.005){
        row = 4
        md.db_05 = c(md.db_05, round(as.double(as.character(df.results[pat,"MD.db"])), 2))
      }
      else if(md == 0.01){
        row = 3
        md.db_1 = c(md.db_1, round(as.double(as.character(df.results[pat,"MD.db"])), 2))
      }
      else if(md == 0.02){
        row = 3
        md.db_1 = c(md.db_1, round(as.double(as.character(df.results[pat,"MD.db"])), 2))
      }
      else if(md == 0.05){
        row = 2
        md.db_5 = c(md.db_5, round(as.double(as.character(df.results[pat,"MD.db"])), 2))
      }
      else if(md == 0.1){
        row = 2
        md.db_5 = c(md.db_5, round(as.double(as.character(df.results[pat,"MD.db"])), 2))
      }
      df.results_md_cluster[row,"Total"] = df.results_md_cluster[row,"Total"] + 1

      for(col in names_col){
        #print(pat)
        if(df.results[pat,col] == T)
          df.results_md_cluster[row,col] = df.results_md_cluster[row,col] + 1

      }
    }
    assign("df.results_md_cluster", df.results_md_cluster, envir = .GlobalEnv)

    # confidence intervals
    CI1 = binom.confint(x=df.results_md_cluster[,2:5], n=df.results_md_cluster[,"Total"], methods="wilson")
    CI2 = binom.confint(x=df.results_md_cluster[,6], n=df.results_md_cluster[,"Total"], methods="wilson")

    # modify name of each bin
    md.db = list(md.db_10, md.db_5, md.db_1, md.db_05)
    for(row in 1:nrow(df.results_md_cluster)){
      df.results_md_cluster[row,1] = paste0(df.results_md_cluster[row,1], "\nmed=", round(median(md.db[[row]]),1), "dB", "\nN=",df.results_md_cluster[row,"Total"])
    }

    # generate data frame to be graphed
    df.results_graph = cbind(MD=df.results_md_cluster[,1], CI1[,7:10], CI2[,4])
    colnames(df.results_graph)[2:6] = names_col[1:5]
    df.results_graph = within(df.results_graph,  MD <- factor(MD, levels=MD))
    print(df.results_graph)

    # melt
    df.melted = melt(df.results_graph, variable.name = "criterion", value.name = "Hit.Rate")
    # round hit rates to 2 sig figs
    df.melted[,"Hit.Rate"] = round(df.melted[,"Hit.Rate"], digits=2)
    # add CI columns
    df.melted = cbind(df.melted,
                      lower.CI=c(CI1[,"lower.HAP2"], CI1[,"lower.UKGTS"], CI1[,"lower.GHT"], CI1[,"lower.FOST"], CI2[,"lower"]),
                      upper.CI=c(CI1[,"upper.HAP2"], CI1[,"upper.UKGTS"], CI1[,"upper.GHT"], CI1[,"upper.FOST"], CI2[,"upper"]))
    print(df.melted)
    plot.hist = ggplot(df.melted, aes(x=MD, y=Hit.Rate, fill=criterion)) +
      geom_linerange(position=position_dodge(0.5), aes(ymin=lower.CI, ymax=upper.CI), size=0.8) +
      geom_point(position = position_dodge(0.5), stat = "identity", aes(fill = criterion), size = 5, shape = 21, colour = "black", size = 5, stroke = 1.3) +
      scale_fill_manual(values = c(GHT = "#bbbcbe", FOST = "#ffffb1", HAP2 = "#ffb1b1", UKGTS = "#b1e6fa", LOGTS = "white")) +
      #geom_text(aes(label = Repeat.Error.Rate, group = criterion), size=6, hjust=0.5, vjust=3, position=position_dodge(0.9)) +
      theme_bw(base_size = 22) +
      ylim(0,1)
    #theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

    #plot.hist = ggplot(df.melted, aes(x=MD, y=Hit.Rate, fill=criterion)) +
    #  geom_bar(position = position_dodge(), stat = "identity") +
    #  scale_fill_manual("criteria", values = c("GHT" = "#bbbcbe", "FOST" = "#ffffb1", "HAP2" = "#ffb1b1", "UKGTS" = "#b1e6fa")) +
    #  geom_errorbar(position=position_dodge(0.9), width=.5, aes(ymin=lower.CI, ymax=upper.CI)) +
    #  geom_text(aes(label = Hit.Rate, group = criterion), size=6, hjust=0.5, vjust=3, position=position_dodge(0.9)) +
    #  theme_bw(base_size = 22) +
    #  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  }
  else if(x_var == "rmd"){

    # confidence intervals
    CI = binom.confint(x=df.criteria_reproducibility_md[,c("HAP2","UKGTS","GHT","FOST","LOGTS")],
                       n=df.criteria_reproducibility_md[,c("HAP2.T","UKGTS.T","GHT.T","FOST.T","LOGTS.T")], methods="wilson")

    # modify name of each bin
    x_names=c("","","","")
    for(row in 1:nrow(df.results_md_cluster)){
      x_names[row] = paste0(df.results_md_cluster[row,1], "\n",
                            df.criteria_reproducibility_md[row,"HAP2.T"], "|",
                            df.criteria_reproducibility_md[row,"UKGTS.T"], "|",
                            df.criteria_reproducibility_md[row,"GHT.T"], "|",
                            df.criteria_reproducibility_md[row,"FOST.T"], "|",
                            df.criteria_reproducibility_md[row,"LOGTS.T"])
    }

    # generate data frame to be graphed
    df.results_graph = cbind(MD=x_names, CI[,12:16])
    colnames(df.results_graph)[2:6] = names_col[1:5]
    df.results_graph = within(df.results_graph,  MD <- factor(MD, levels=MD))
    print(df.results_graph)

    # melt
    df.melted = melt(df.results_graph, variable.name = "criterion", value.name = "Repeat.Error.Rate")
    # round hit rates to 2 sig figs
    df.melted[,"Repeat.Error.Rate"] = round(df.melted[,"Repeat.Error.Rate"], digits=2)
    # add CI columns
    df.melted = cbind(df.melted,
                      lower.CI=c(CI[,"lower.HAP2"], CI[,"lower.UKGTS"], CI[,"lower.GHT"], CI[,"lower.FOST"], CI[,"lower.LOGTS"]),
                      upper.CI=c(CI[,"upper.HAP2"], CI[,"upper.UKGTS"], CI[,"upper.GHT"], CI[,"upper.FOST"], CI[,"upper.LOGTS"]))
    print(df.melted)
    plot.hist = ggplot(df.melted, aes(x=MD, y=Repeat.Error.Rate, fill=criterion)) +
      geom_linerange(position=position_dodge(0.5), aes(ymin=lower.CI, ymax=upper.CI), size=0.8) +
      geom_point(position = position_dodge(0.5), stat = "identity", aes(fill = criterion), size = 5, shape = 21, colour = "black", size = 5, stroke = 1.3) +
      scale_fill_manual(values = c(GHT = "#bbbcbe", FOST = "#ffffb1", HAP2 = "#ffb1b1", UKGTS = "#b1e6fa", LOGTS = "white")) +
      #geom_text(aes(label = Repeat.Error.Rate, group = criterion), size=6, hjust=0.5, vjust=3, position=position_dodge(0.9)) +
      theme_bw(base_size = 22) +
      ylim(0,1)
  }
  else if(x_var == "ros"){

    # confidence intervals
    CI = binom.confint(x=df.criteria_reproducibility_os[,c("HAP2","UKGTS","GHT","FOST","LOGTS")],
                        n=df.criteria_reproducibility_os[,c("HAP2.T","UKGTS.T","GHT.T","FOST.T","LOGTS.T")], methods="wilson")

    x_names=c("","","","")
    for(row in 1:nrow(df.results_md_cluster)){
      x_names[row] = paste0(df.results_os_cluster[row,1], "\n",
                            df.criteria_reproducibility_os[row,"HAP2.T"], "|",
                            df.criteria_reproducibility_os[row,"UKGTS.T"], "|",
                            df.criteria_reproducibility_os[row,"GHT.T"], "|",
                            df.criteria_reproducibility_os[row,"FOST.T"], "|",
                            df.criteria_reproducibility_os[row,"LOGTS.T"])
    }

    # generate data frame to be graphed
    df.results_graph = cbind(OCT.Score=x_names, CI[,12:16])
    colnames(df.results_graph)[2:6] = names_col[1:5]
    df.results_graph = within(df.results_graph,  OCT.Score <- factor(OCT.Score, levels=OCT.Score))
    print(df.results_graph)

    # melt
    df.melted = melt(df.results_graph, variable.name = "criterion", value.name = "Repeat.Error.Rate")
    # round hit rates to 2 sig figs
    df.melted[,"Repeat.Error.Rate"] = round(df.melted[,"Repeat.Error.Rate"], digits=2)
    # add CI columns
    df.melted = cbind(df.melted,
                      lower.CI=c(CI[,"lower.HAP2"], CI[,"lower.UKGTS"], CI[,"lower.GHT"], CI[,"lower.FOST"], CI[,"lower.LOGTS"]),
                      upper.CI=c(CI[,"upper.HAP2"], CI[,"upper.UKGTS"], CI[,"upper.GHT"], CI[,"upper.FOST"], CI[,"upper.LOGTS"]))
    print(df.melted)
    plot.hist = ggplot(df.melted, aes(x=OCT.Score, y=Repeat.Error.Rate, fill=criterion)) +
      geom_linerange(position=position_dodge(0.5), aes(ymin=lower.CI, ymax=upper.CI), size=0.8) +
      geom_point(position = position_dodge(0.5), stat = "identity", aes(fill = criterion), size = 5, shape = 21, colour = "black", size = 5, stroke = 1.3) +
      scale_fill_manual(values = c(GHT = "#bbbcbe", FOST = "#ffffb1", HAP2 = "#ffb1b1", UKGTS = "#b1e6fa", LOGTS = "white")) +
      #geom_text(aes(label = Repeat.Error.Rate, group = criterion), size=6, hjust=0.5, vjust=3, position=position_dodge(0.9)) +
      theme_bw(base_size = 22) +
      ylim(0,1)

  }

  return(plot.hist)
}

#' Euler & Clustered Histo for MD subgroups
#'
#' Print Euler diagrams and Clustered Histogram for each fo the MD percentiles (10, 5, 2, 1, and 0.5 %)
#' @export
printFigure1 = function(df.results=df.best_match)
{
  names = c("FOST", "HAP2", "UKGTS", "GHT", "LOGTS",
            "FOST&HAP2", "FOST&UKGTS", "FOST&GHT", "FOST&LOGTS", "HAP2&UKGTS",
            "HAP2&GHT", "HAP2&LOGTS", "UKGTS&GHT", "UKGTS&LOGTS", "GHT&LOGTS",
            "FOST&HAP2&UKGTS", "FOST&HAP2&GHT", "FOST&HAP2&LOGTS", "FOST&UKGTS&GHT", "FOST&UKGTS&LOGTS",
            "FOST&GHT&LOGTS", "HAP2&UKGTS&GHT", "HAP2&GHT&LOGTS", "HAP2&UKGTS&LOGTS", "UKGTS&GHT&LOGTS",
            "FOST&HAP2&UKGTS&GHT", "FOST&HAP2&UKGTS&LOGTS", "FOST&HAP2&GHT&LOGTS", "FOST&UKGTS&GHT&LOGTS", "HAP2&UKGTS&GHT&LOGTS",
            "FOST&HAP2&UKGTS&GHT&LOGTS", "Neither")

  md_subgroups = c("<0.5%", "0.5-2%", "2-10%", ">10.0%")
  df.results_md_05 = df.results_md_2 = df.results_md_10 = df.results_md_11 = data.frame(FOST=logical(), GHT=logical(), HAP2=logical(), LOGTS=logical(), UKGTS=logical())
  df.interx = data.frame(matrix(NA, nrow = length(names), ncol = length(md_subgroups)))
  rownames(df.interx) = names
  colnames(df.interx) = md_subgroups
  assign("df.interx", df.interx, envir = .GlobalEnv)


  col_range = which(colnames(df.results)=="GHT"):which(colnames(df.results)=="UKGTS")
  for(pat in 1:nrow(df.results)){
    md = df.results[pat,"MD.pval"]
    if(is.na(md)){
      df.results_md_11 = rbind(df.results_md_11, df.results[pat,col_range])
    }
    else if(md == 0.005){
      df.results_md_05 = rbind(df.results_md_05, df.results[pat,col_range])
    }
    else if(md == 0.01){
      df.results_md_2 = rbind(df.results_md_2, df.results[pat,col_range])
    }
    else if(md == 0.02){
      df.results_md_2 = rbind(df.results_md_2, df.results[pat,col_range])
    }
    else if(md == 0.05){
      df.results_md_10 = rbind(df.results_md_10, df.results[pat,col_range])
    }
    else if(md == 0.1){
      df.results_md_10 = rbind(df.results_md_10, df.results[pat,col_range])
    }
  }
  ls.results_md = list(df.results_md_05, df.results_md_2, df.results_md_10, df.results_md_11)

  par(mfrow=c(2,2), mar=c(0.2,0.2,0.2,0.2))
  for(i in 1:length(ls.results_md)){
    printEulerDiag(md_subgroups[i], as.data.frame(ls.results_md[i]))
  }

  p = recordPlot()
  dev.off()
  h = printClusteredHist(df.results, "md")
  plot_grid(p, h, ncol = 1, labels = c("", ""))
}

#' Euler & Clustered Histo for OCT Score subgroups
#'
#' Print Euler diagrams and Clustered Histogram for each fo the OCT scores (0-6)
#' @export
printFigure2 = function(df.results=df.best_match)
{
  names = c("FOST", "HAP2", "UKGTS", "GHT", "LOGTS",
            "FOST&HAP2", "FOST&UKGTS", "FOST&GHT", "FOST&LOGTS", "HAP2&UKGTS",
            "HAP2&GHT", "HAP2&LOGTS", "UKGTS&GHT", "UKGTS&LOGTS", "GHT&LOGTS",
            "FOST&HAP2&UKGTS", "FOST&HAP2&GHT", "FOST&HAP2&LOGTS", "FOST&UKGTS&GHT", "FOST&UKGTS&LOGTS",
            "FOST&GHT&LOGTS", "HAP2&UKGTS&GHT", "HAP2&GHT&LOGTS", "HAP2&UKGTS&LOGTS", "UKGTS&GHT&LOGTS",
            "FOST&HAP2&UKGTS&GHT", "FOST&HAP2&UKGTS&LOGTS", "FOST&HAP2&GHT&LOGTS", "FOST&UKGTS&GHT&LOGTS", "HAP2&UKGTS&GHT&LOGTS",
            "FOST&HAP2&UKGTS&GHT&LOGTS", "Neither")

  os_subgroups = c("6", "4-5", "1-3", "0")
  df.results_os_0 = df.results_os_13 = df.results_os_45 = df.results_os_6 = data.frame(FOST=logical(), GHT=logical(), HAP2=logical(), LOGTS=logical(), UKGTS=logical())
  df.interx = data.frame(matrix(NA, nrow = length(names), ncol = length(os_subgroups)))
  rownames(df.interx) = names
  colnames(df.interx) = os_subgroups
  assign("df.interx", df.interx, envir = .GlobalEnv)


  col_range = which(colnames(df.results)=="GHT"):which(colnames(df.results)=="UKGTS")
  for(pat in 1:nrow(df.results)){
    os = df.results[pat,"OCT.Score"]
    if(os == 0){
      df.results_os_0 = rbind(df.results_os_0, df.results[pat,col_range])
    }
    else if(inside.range(os, c(1,3))){
      df.results_os_13 = rbind(df.results_os_13, df.results[pat,col_range])
    }
    else if(inside.range(os, c(4,5))){
      df.results_os_45 = rbind(df.results_os_45, df.results[pat,col_range])
    }
    else if(os == 6){
      df.results_os_6 = rbind(df.results_os_6, df.results[pat,col_range])
    }
  }
  ls.results_os = list(df.results_os_6, df.results_os_45, df.results_os_13, df.results_os_0)

  par(mfrow=c(2,2), mar=c(0.2,0.2,0.2,0.2))
  for(i in 1:length(ls.results_os)){
    printEulerDiag(os_subgroups[i], as.data.frame(ls.results_os[i]))
  }

  p = recordPlot()
  dev.off()
  h = printClusteredHist(df.results, "os")
  plot_grid(p, h, ncol = 1, labels = c("", ""))
}

#' Print to pdf
#'
#' Create a .pdf file with desired STATPAC SFA plots.
#' @export
createPdf = function(first_pat = 1, last_pat = NUM_PAT)
{
   if((first_pat < 1) || (first_pat > NUM_PAT))
    print("Invaild first patient number.")
  else if((last_pat < 1) || (last_pat > NUM_PAT))
    print("Invaild last patient number.")
  else if(first_pat > last_pat)
    print("Last patient number must be larger than first patient number.")
  else {
    pdf("test_vf_maps.pdf")

    for (pat in first_pat:last_pat){
      plot(printVfMap(pat))
    }
    dev.off()
  }
}


#' VF criteria application
#'
#' Apply Foster, HAP2, UKGTS criteria to all patients from a properly-formatted excel spreadsheet
#' @export
assignVfCriteria = function()
{
  assign("df.criteria_results", data.frame("Patient.ID"=integer(), "Eye"=integer(), "Date.Time"=integer(),
                                           "MD.db"=double(), "MD.pval"=double(),
                                           "GHT"=integer(), "FOST"=integer(), "HAP2"=integer(),
                                           "LOGTS"=integer(), "UKGTS"=integer(), "EAGLE"=integer(), "AGIS"=integer()),
                                           envir = .GlobalEnv)
  for(pat in 1:NUM_PAT)
  {
    df.criteria_results[pat,] <<- c(df.vf_data[(pat+VF_V_OFST),1], df.vf_data[(pat+VF_V_OFST),4],
                                    df.vf_data[(pat+VF_V_OFST),3], df.vf_data[(pat+VF_V_OFST),MD_H_OFST-1],
                                    df.vf_data[(pat+VF_V_OFST),MD_H_OFST],
                                    checkGhtCriteria(pat), checkFostCriteria(pat), checkMhpaCriteria(pat),
                                    checkLogtsCriteria(pat), checkUkgtsCriteria(pat), checkEagleCriteria(pat), checkAgisCriteria(pat))
  }
}

#' GHT criteria (GHT)
#'
#' Glaucoma Hemifield Test (GHT) “outside normal limits”
#' @param pat_id Patient position in the list on the excel spreadsheet
#' @return BOOL Whether patient VF analysis resulted in a positive result.
#' @export
checkGhtCriteria = function(pat_id)
{
  ght = df.vf_data[VF_V_OFST+pat_id,GHT_H_OFST]
  if(ght == "Outside normal limits")
    return(TRUE)

  return(FALSE)
}

#' Foster criteria (FOST)
#'
#' Glaucoma Hemifield Test (GHT) “outside normal limits” AND a cluster of three contiguous points at the 5% level on the pattern deviation plot
#' @param pat_id Patient position in the list on the excel spreadsheet
#' @return BOOL Whether patient VF analysis resulted in a positive result.
#' @export
checkFostCriteria = function(pat_id)
{
  ght = df.vf_data[VF_V_OFST+pat_id,GHT_H_OFST]
  if(ght != "Outside normal limits")
    return(FALSE)

  df.VF_PDp = data.frame(t(df.vf_data[VF_V_OFST+pat_id, PDP_H_OFST:(PDP_H_OFST+VF_H_SIZE-1)]))


  for(pt1 in 1:VF_H_SIZE)
  {
    if(!(is.na(df.VF_PDp[pt1,1])))
    {
      if((df.VF_PDp[pt1,1] <= 5))
      {
        print(paste0("pt1=", pt1))

        for(n1 in 4:11)
        {
          pt2 = df.pt_mapping[pt1,n1]

          if(!(is.na(df.VF_PDp[pt2,1])))
          {
            if(df.VF_PDp[pt2,1] <= 5)
            {
              print(paste0("pt2=", pt2))

              for(n2 in 4:11)
              {
                pt3 = df.pt_mapping[pt2,n2]

                if(!(is.na(df.VF_PDp[pt3,1])))
                {
                  if(df.VF_PDp[pt3,1] <= 5)
                  {
                    if(pt1 != pt3)
                    {
                      print(paste0("pt3=", pt3))
                      return(TRUE)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(FALSE)
}

#' Modified Hoddap-Parrish-Anderson criteria (HAP2)
#'
#' GHT outside normal limits; OR cluster of 3 points on the pattern deviation plot depressed at p<5%, one of which depressed at p<1%; OR PSD with p<5%
#' @param pat_id Patient position in the list on the excel spreadsheet
#' @return BOOL Whether patient VF analysis resulted in a positive result.
#' @export
checkMhpaCriteria = function(pat_id)
{
  # check GHT
  ght = df.vf_data[VF_V_OFST+pat_id,GHT_H_OFST]
  if(ght == "Outside normal limits")
    return(TRUE)

  # check PSD
  psd = df.vf_data[VF_V_OFST+pat_id,PSD_H_OFST]
  if(!(is.na(psd)))
  {
    if(psd < 0.1)
      return(TRUE)
  }

  # check cluster of 3 points
  df.VF_PDp = data.frame(t(df.vf_data[VF_V_OFST+pat_id, PDP_H_OFST:(PDP_H_OFST+VF_H_SIZE-1)]))

  for(pt1 in 1:VF_H_SIZE)
  {
    if(!(is.na(df.VF_PDp[pt1,1])))
    {
      if((df.VF_PDp[pt1,1] <= 5))
      {
        print(paste0("pt1=", pt1))

        for(n1 in 4:11)
        {
          pt2 = df.pt_mapping[pt1,n1]

          if(!(is.na(df.VF_PDp[pt2,1])))
          {
            if(df.VF_PDp[pt2,1] <= 5)
            {
              print(paste0("pt2=", pt2))

              for(n2 in 4:11)
              {
                pt3 = df.pt_mapping[pt2,n2]

                if(!(is.na(df.VF_PDp[df.pt_mapping[pt2,n2],1])))
                {
                  if(df.VF_PDp[df.pt_mapping[pt2,n2],1] <= 5)
                  {
                    if(pt1 != pt3)
                    {
                      print(paste0("pt3=", pt3))

                      if((df.VF_PDp[pt1,1] <= 1) || (df.VF_PDp[pt2,1] <= 1) || (df.VF_PDp[pt3,1] <= 1))
                        return(TRUE)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(FALSE)
}

#' UKGTS criteria
#'
#' 2 or more contiguous points with P<0.01 loss or more; OR 3 or more contiguous points with P<0.05 loss or more; OR a 10-dB difference across the nasal horizontal midline at 2 or more adjacent points in the total deviation plot.
#' @param pat_id Patient position in the list on the excel spreadsheet
#' @return BOOL Whether patient VF analysis resulted in a positive result.
#' @export
checkUkgtsCriteria = function(pat_id)
{
  #import data points
  df.VF_TD = data.frame(t(df.vf_data[VF_V_OFST+pat_id, TD_H_OFST:(TD_H_OFST+VF_H_SIZE-1)]))
  df.VF_TDp = data.frame(t(df.vf_data[VF_V_OFST+pat_id, TDP_H_OFST:(TDP_H_OFST+VF_H_SIZE-1)]))

  # check cluster of 2 points
  for(pt1 in 1:VF_H_SIZE)
  {
    if(!(is.na(df.VF_TDp[pt1,1])))
    {
      if((df.VF_TDp[pt1,1] <= 1))
      {
        print(paste0("pt1=", pt1))

        for(n1 in 4:11)
        {
          pt2 = df.pt_mapping[pt1,n1]

          if(!(is.na(df.VF_TDp[pt2,1])))
          {
            if(df.VF_TDp[pt2,1] <= 1)
            {
              print(paste0("pt2=", pt2))

              return(TRUE)
            }
          }
        }
      }
    }
  }

  # check cluster of 3 points
  for(pt1 in 1:VF_H_SIZE)
  {
    if(!(is.na(df.VF_TDp[pt1,1])))
    {
      if((df.VF_TDp[pt1,1] <= 5))
      {
        print(paste0("pt1=", pt1))

        for(n1 in 4:11)
        {
          pt2 = df.pt_mapping[pt1,n1]

          if(!(is.na(df.VF_TDp[pt2,1])))
          {
            if(df.VF_TDp[pt2,1] <= 5)
            {
              print(paste0("pt2=", pt2))

              for(n2 in 4:11)
              {
                pt3 = df.pt_mapping[pt2,n2]

                if(!(is.na(df.VF_TDp[pt3,1])))
                {
                  if(df.VF_TDp[pt3,1] <= 5)
                  {
                    if(pt1 != pt3)
                    {
                      print(paste0("pt3=", pt3))
                      return(TRUE)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  # check reflection points
  for(pt1 in 1:(VF_H_SIZE/2))
  {
    pt1_r = df.pt_mapping[pt1,12]

    if((!(is.na(df.VF_TD[pt1,1]))) && (!(is.na(df.VF_TD[pt1_r,1]))))
    {
      if(abs(df.VF_TD[pt1,1] - df.VF_TD[pt1_r,1]) >= 10)
      {
        print(paste0("pt1=", pt1, " pt1_r=", pt1_r))

        for(n1 in 4:11)
        {
          pt2 = df.pt_mapping[pt1,n1]
          pt2_r = df.pt_mapping[pt2,12]

          if((!(is.na(df.VF_TD[pt2,1]))) && (!(is.na(df.VF_TD[pt2_r,1]))))
          {
            if(abs(df.VF_TD[pt2,1] - df.VF_TD[pt2_r,1]) >= 10)
            {
              if((pt1 != pt2) && (pt1_r != pt2))
              {
                print(paste0("pt2=", pt2, " pt2_r=", pt2_r))
                return(TRUE)
              }

            }
          }
        }
      }
    }
  }
  return(FALSE)
}

#' Low-pressure glaucoma treatment study criteria (LOGTS)
#'
#' Presence of at least 3 contiguous points depressed more than 8 decibels or 2 contiguous points depressed more than 10 decibels on TD plot
#' @param pat_id Patient position in the list on the excel spreadsheet
#' @return BOOL Whether patient VF analysis resulted in a positive result.
#' @export
checkLogtsCriteria = function(pat_id)
{
  #import data points
  df.VF_TD = data.frame(t(df.vf_data[VF_V_OFST+pat_id, TD_H_OFST:(TD_H_OFST+VF_H_SIZE-1)]))

  # check cluster of 2 points
  for(pt1 in 1:VF_H_SIZE)
  {
    if(!(is.na(df.VF_TD[pt1,1])))
    {
      if((df.VF_TD[pt1,1] < -10))
      {
        print(paste0("pt1=", pt1))

        for(n1 in 4:11)
        {
          pt2 = df.pt_mapping[pt1,n1]

          if(!(is.na(df.VF_TD[pt2,1])))
          {
            if(df.VF_TD[pt2,1] < -10)
            {
              print(paste0("pt2=", pt2))

              return(TRUE)
            }
          }
        }
      }
    }
  }

  # check cluster of 3 points
  for(pt1 in 1:VF_H_SIZE)
  {
    if(!(is.na(df.VF_TD[pt1,1])))
    {
      if((df.VF_TD[pt1,1] < -8))
      {
        print(paste0("pt1=", pt1))

        for(n1 in 4:11)
        {
          pt2 = df.pt_mapping[pt1,n1]

          if(!(is.na(df.VF_TD[pt2,1])))
          {
            if(df.VF_TD[pt2,1] < -8)
            {
              print(paste0("pt2=", pt2))

              for(n2 in 4:11)
              {
                pt3 = df.pt_mapping[pt2,n2]

                if(!(is.na(df.VF_TD[pt3,1])))
                {
                  if(df.VF_TD[pt3,1] < -8)
                  {
                    if(pt1 != pt3)
                    {
                      print(paste0("pt3=", pt3))
                      return(TRUE)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(FALSE)
}

checkEagleCriteria = function(pat_id)
{
  return(FALSE)
}

checkAgisCriteria = function(pat_id)
{
  return(FALSE)
}

#' Match VF to OCT
#'
#' Matches VF criteria results to OCT data, selecting any first exam above OCT quality of 15
#' @param window Maximum number of months time difference between VF and OCT exams
#' @return none
#' @export
matchVfToOct = function(window = 4)
{
  assign("df.criteria_results", arrange(df.criteria_results, Patient.ID), envir = .GlobalEnv)
  assign("df.rnfl_data", arrange(df.rnfl_data, StudyID), envir = .GlobalEnv)
  assign("df.mrw_data", arrange(df.mrw_data, Patient.ID), envir = .GlobalEnv)


  assign("df.gcl_data", arrange(df.gcl_data, StudyID), envir = .GlobalEnv)

  ls.criteria_results.unique = unique(df.criteria_results[,"Patient.ID"])

  names = c("Patient.ID",
            "Eye.Vf", "Eye.RNFL", "Eye.MRW", "Eye.GCL",
            "Date.Vf", "Date.RNFL", "Date.MRW", "Date.GCL",
            "Quality.RNFL", "Quality.MRW", "Quality.GCL",
            "GHT", "FOST", "HAP2", "LOGTS", "UKGTS", "EAGLE", "AGIS",
            "MD.db", "MD.pval",
            "RNFLClass_G", "RNFLClass_T", "RNFLClass_TS", "RNFLClass_TI", "RNFLClass_N", "RNFLClass_NS", "RNFLClass_NI",
            "MRW.P.Global", "MRW.P.Tmp", "MRW.P.TS", "MRW.P.TI", "MRW.P.Nas", "MRW.P.NS", "MRW.P.NI",
            "GCL.ABN.Global", "GCL.ABN.Sup", "GCL.ABN.TS", "GCL.ABN.NS", "GCL.ABN.Inf", "GCL.ABN.NI", "GCL.ABN.TI")

  df.best_match = data.frame(matrix(NA, nrow = 0, ncol = length(names)))
  colnames(df.best_match) = names

  time_bound = 2592000*window

  for(vf.v in ls.criteria_results.unique){
    mtx.rnfl = matrix(0, nrow = 1, ncol = 4)
    colnames(mtx.rnfl) = c("vf.i", "i", "eye", "q")
    mtx.mrw = matrix(0, nrow = 1, ncol = 4)
    colnames(mtx.mrw) = c("vf.i", "i", "eye", "q")
    mtx.gcl = matrix(0, nrow = 1, ncol = 4)
    colnames(mtx.gcl) = c("vf.i", "i", "eye", "q")

    mtx.max = matrix(integer(), nrow = 0, ncol = 7)
    colnames(mtx.max) = c("vf.i", "rnfl.i", "mrw.i", "gcl.i","rnfl.q", "mrw.q", "gcl.q")

    for(vf.i in which(df.criteria_results[,"Patient.ID"] == vf.v)){
      mtx.rnfl[1,] = c(0,0,0,0)
      mtx.mrw[1,]  = c(0,0,0,0)
      mtx.gcl[1,]  = c(0,0,0,0)

      if(df.criteria_results[vf.i,"Eye"] == "L")
        vf.eye = 1
      else
        vf.eye = 2
      vf.epoch = as.integer(df.criteria_results[vf.i,"Date.Time"])
      vf.epoch.upper_bound = vf.epoch + time_bound
      vf.epoch.lower_bound = vf.epoch - time_bound

      for(rnfl.i in which(df.rnfl_data[,"StudyID"] == vf.v)){
        rnfl.eye = as.integer(df.rnfl_data[rnfl.i,"Eye"])
        rnfl.epoch = as.integer(as.POSIXct(df.rnfl_data[rnfl.i,"ExamDate"]))

        if((rnfl.eye == vf.eye) && inside.range(rnfl.epoch, c(vf.epoch.lower_bound, vf.epoch.upper_bound))){
          mtx.rnfl[1,] = c(vf.i, rnfl.i, rnfl.eye, df.rnfl_data[rnfl.i,"Quality"])
          break
        }
      }

      for(mrw.i in which(df.mrw_data[,"Patient.ID"] == vf.v)){
        mrw.eye = as.integer(df.mrw_data[mrw.i,"Eye"])
        mrw.epoch = as.integer(as.POSIXct(df.mrw_data[mrw.i,"ExamDate"]))

        if((mrw.eye == vf.eye) && inside.range(mrw.epoch, c(vf.epoch.lower_bound, vf.epoch.upper_bound))){
          mtx.mrw[1,] = c(vf.i, mrw.i, mrw.eye, df.mrw_data[mrw.i,"Mean.Quality"])
          break
        }
      }

      for(gcl.i in which(df.gcl_data[,"StudyID"] == vf.v)){
        gcl.eye = as.integer(df.gcl_data[gcl.i,"Eye"])
        gcl.epoch = df.gcl_data[gcl.i,"ExamDate"]

        if((gcl.eye == vf.eye) && inside.range(gcl.epoch, c(vf.epoch.lower_bound, vf.epoch.upper_bound))){
          mtx.gcl[1,] = c(vf.i, gcl.i, gcl.eye, df.gcl_data[gcl.i,"Quality"])
          break
        }
      }
      if((mtx.rnfl[1,1] != 0) & (mtx.mrw[1,1] != 0) & (mtx.gcl[1,1] != 0)){
        break
      }
    }

    if((mtx.rnfl[1,1] != 0) & (mtx.mrw[1,1] != 0) & (mtx.gcl[1,1] != 0)){
      vf.idx   = mtx.rnfl[1,"vf.i"]
      rnfl.idx = mtx.rnfl[1,"i"]
      mrw.idx  = mtx.mrw[1,"i"]
      gcl.idx  = mtx.gcl[1,"i"]

      best_match = data.frame(df.criteria_results[vf.idx,"Patient.ID"],
                              df.criteria_results[vf.idx, "Eye"], df.rnfl_data[rnfl.idx,"Eye"], df.mrw_data[mrw.idx,"Eye"], df.gcl_data[gcl.idx,"Eye"],
                              as.Date.POSIXct(as.integer(df.criteria_results[vf.idx,"Date.Time"])),
                              df.rnfl_data[rnfl.idx,"ExamDate"], df.mrw_data[mrw.idx,"ExamDate"], as.Date.POSIXct(df.gcl_data[gcl.idx,"ExamDate"]),
                              df.rnfl_data[rnfl.idx,"Quality"], df.mrw_data[mrw.idx,"Mean.Quality"], df.gcl_data[gcl.idx,"Quality"],
                              df.criteria_results[vf.idx, 6:12], df.criteria_results[vf.idx,"MD.db"], df.criteria_results[vf.idx,"MD.pval"],
                              df.rnfl_data[rnfl.idx, 19:25],
                              df.mrw_data[mrw.idx, 134:140],
                              df.gcl_data[gcl.idx, 4:10])

      colnames(best_match) = names

      df.best_match = rbind(df.best_match, best_match)
    }
  }

  assign("df.best_match", df.best_match, envir = .GlobalEnv)
}

#' Match VF to OCT 2
#'
#' Matches VF criteria results to OCT data, selecting the highest quality OCT exams
#' @param window Maximum number of months time difference between VF and OCT exams
#' @return none
#' @export
matchVfToOct2 = function(window = 4)
{
  assign("df.criteria_results", arrange(df.criteria_results, Patient.ID), envir = .GlobalEnv)
  assign("df.rnfl_data", arrange(df.rnfl_data, StudyID), envir = .GlobalEnv)
  assign("df.mrw_data", arrange(df.mrw_data, StudyID), envir = .GlobalEnv)


  assign("df.gcl_data", arrange(df.gcl_data, StudyID), envir = .GlobalEnv)

  ls.criteria_results.unique = unique(df.criteria_results[,"Patient.ID"])

  names = c("Patient.ID",
            "Eye.Vf", "Eye.RNFL", "Eye.MRW", "Eye.GCL",
            "Date.Vf", "Date.RNFL", "Date.MRW", "Date.GCL",
            "Quality.RNFL", "Quality.MRW", "Quality.GCL",
            "GHT", "FOST", "HAP2", "LOGTS", "UKGTS", "EAGLE", "AGIS",
            "MD.db", "MD.pval",
            "RNFLClass_G", "RNFLClass_T", "RNFLClass_TS", "RNFLClass_TI", "RNFLClass_N", "RNFLClass_NS", "RNFLClass_NI",
            "MRW.P.Global", "MRW.P.Tmp", "MRW.P.TS", "MRW.P.TI", "MRW.P.Nas", "MRW.P.NS", "MRW.P.NI",
            "GCL.ABN.Global", "GCL.ABN.Sup", "GCL.ABN.TS", "GCL.ABN.NS", "GCL.ABN.Inf", "GCL.ABN.NI", "GCL.ABN.TI")

  df.best_match = data.frame(matrix(NA, nrow = 0, ncol = length(names)))
  colnames(df.best_match) = names

  time_bound = 2592000*window

  for(vf.v in ls.criteria_results.unique){
    mtx.rnfl = matrix(integer(), nrow = 0, ncol = 4)
    colnames(mtx.rnfl) = c("vf.i", "i", "eye", "q")
    mtx.mrw = matrix(integer(), nrow = 0, ncol = 4)
    colnames(mtx.mrw) = c("vf.i", "i", "eye", "q")
    mtx.gcl = matrix(integer(), nrow = 0, ncol = 4)
    colnames(mtx.gcl) = c("vf.i", "i", "eye", "q")

    mtx.max = matrix(integer(), nrow = 0, ncol = 7)
    colnames(mtx.max) = c("vf.i", "rnfl.i", "mrw.i", "gcl.i","rnfl.q", "mrw.q", "gcl.q")

    for(vf.i in which(df.criteria_results[,"Patient.ID"] == vf.v)){
      if(df.criteria_results[vf.i,"Eye"] == "L")
        vf.eye = 1
      else
        vf.eye = 2
      vf.epoch = as.integer(df.criteria_results[vf.i,"Date.Time"])
      vf.epoch.upper_bound = vf.epoch + time_bound
      vf.epoch.lower_bound = vf.epoch - time_bound

      for(rnfl.i in which(df.rnfl_data[,"StudyID"] == vf.v)){
        rnfl.eye = as.integer(df.rnfl_data[rnfl.i,"Eye"])
        rnfl.epoch = as.integer(as.POSIXct(df.rnfl_data[rnfl.i,"ExamDate"]))

        if((rnfl.eye == vf.eye) && inside.range(rnfl.epoch, c(vf.epoch.lower_bound, vf.epoch.upper_bound))){
          mtx.rnfl = rbind(mtx.rnfl, c(vf.i, rnfl.i, rnfl.eye, df.rnfl_data[rnfl.i,"Quality"]))
        }
      }

      for(mrw.i in which(df.mrw_data[,"StudyID"] == vf.v)){
        mrw.eye = as.integer(df.mrw_data[mrw.i,"Eye"])
        mrw.epoch = as.integer(as.POSIXct(df.mrw_data[mrw.i,"ExamDate"]))

        if((mrw.eye == vf.eye) && inside.range(mrw.epoch, c(vf.epoch.lower_bound, vf.epoch.upper_bound))){
          mtx.mrw = rbind(mtx.mrw, c(vf.i, mrw.i, mrw.eye, df.mrw_data[mrw.i,"Mean.Quality"]))
        }
      }

      for(gcl.i in which(df.gcl_data[,"StudyID"] == vf.v)){
        gcl.eye = as.integer(df.gcl_data[gcl.i,"Eye"])
        gcl.epoch = df.gcl_data[gcl.i,"ExamDate"]

        if((gcl.eye == vf.eye) && inside.range(gcl.epoch, c(vf.epoch.lower_bound, vf.epoch.upper_bound))){
          mtx.gcl = rbind(mtx.gcl, c(vf.i, gcl.i, gcl.eye, df.gcl_data[gcl.i,"Quality"]))
        }
      }
    }


    for(vf.i in unique(mtx.rnfl[,"vf.i"])){
      #if(length(which(mtx.mrw[,"vf.i"] == vf.i)) && length(which(mtx.gcl[,"vf.i"] == vf.i))){
      rnfl.j = which((mtx.rnfl[,"vf.i"] == vf.i) & (mtx.rnfl[,"eye"] == 1))[which.max(mtx.rnfl[which((mtx.rnfl[,"vf.i"] == vf.i) & (mtx.rnfl[,"eye"] == 1)), "q"])]
      mrw.j  = which((mtx.mrw [,"vf.i"] == vf.i) & (mtx.mrw [,"eye"] == 1))[which.max(mtx.mrw [which((mtx.mrw [,"vf.i"] == vf.i) & (mtx.mrw [,"eye"] == 1)), "q"])]
      gcl.j  = which((mtx.gcl [,"vf.i"] == vf.i) & (mtx.gcl [,"eye"] == 1))[which.max(mtx.gcl [which((mtx.gcl [,"vf.i"] == vf.i) & (mtx.gcl [,"eye"] == 1)), "q"])]

      if(length(rnfl.j) && length(mrw.j) && length(gcl.j))
        mtx.max = rbind(mtx.max, c(vf.i, mtx.rnfl[rnfl.j,"i"], mtx.mrw[mrw.j,"i"], mtx.gcl[gcl.j,"i"],
                                   mtx.rnfl[rnfl.j,"q"], mtx.mrw[mrw.j,"q"], mtx.gcl[gcl.j,"q"]))

      rnfl.j = which((mtx.rnfl[,"vf.i"] == vf.i) & (mtx.rnfl[,"eye"] == 2))[which.max(mtx.rnfl[which((mtx.rnfl[,"vf.i"] == vf.i) & (mtx.rnfl[,"eye"] == 2)), "q"])]
      mrw.j  = which((mtx.mrw [,"vf.i"] == vf.i) & (mtx.mrw [,"eye"] == 2))[which.max(mtx.mrw [which((mtx.mrw [,"vf.i"] == vf.i) & (mtx.mrw [,"eye"] == 2)), "q"])]
      gcl.j  = which((mtx.gcl [,"vf.i"] == vf.i) & (mtx.gcl [,"eye"] == 2))[which.max(mtx.gcl [which((mtx.gcl [,"vf.i"] == vf.i) & (mtx.gcl [,"eye"] == 2)), "q"])]

      if(length(rnfl.j) && length(mrw.j) && length(gcl.j))
        mtx.max = rbind(mtx.max, c(vf.i, mtx.rnfl[rnfl.j,"i"], mtx.mrw[mrw.j,"i"], mtx.gcl[gcl.j,"i"],
                                   mtx.rnfl[rnfl.j,"q"], mtx.mrw[mrw.j,"q"], mtx.gcl[gcl.j,"q"]))
      #}
    }
    idx = which.max(as.integer(mtx.max[,"rnfl.q"]) + as.integer(mtx.max[,"mrw.q"]) + as.integer(mtx.max[,"gcl.q"]))
    vf.idx   = mtx.max[idx,"vf.i"]
    rnfl.idx = mtx.max[idx,"rnfl.i"]
    mrw.idx  = mtx.max[idx,"mrw.i"]
    gcl.idx  = mtx.max[idx,"gcl.i"]

    best_match = data.frame(df.criteria_results[vf.idx,"Patient.ID"],
                            df.criteria_results[vf.idx, "Eye"], df.rnfl_data[rnfl.idx,"Eye"], df.mrw_data[mrw.idx,"Eye"], df.gcl_data[gcl.idx,"Eye"],
                            as.Date.POSIXct(as.integer(df.criteria_results[vf.idx,"Date.Time"])),
                            df.rnfl_data[rnfl.idx,"ExamDate"], df.mrw_data[mrw.idx,"ExamDate"], as.Date.POSIXct(df.gcl_data[gcl.idx,"ExamDate"]),
                            df.rnfl_data[rnfl.idx,"Quality"], df.mrw_data[mrw.idx,"Mean.Quality"], df.gcl_data[gcl.idx,"Quality"],
                            df.criteria_results[vf.idx, 6:12], df.criteria_results[vf.idx,"MD.db"], df.criteria_results[vf.idx,"MD.pval"],
                            df.rnfl_data[rnfl.idx, 19:25],
                            df.mrw_data[mrw.idx, 134:140],
                            df.gcl_data[gcl.idx, 4:10])

    colnames(best_match) = names

    df.best_match = rbind(df.best_match, best_match)

    #if(a == 4){
    #  a = 0
    #  break;
    #}
    #else
    #  a = a + 1
  }

  assign("df.best_match", df.best_match, envir = .GlobalEnv)
}

#' OCT Score
#'
#' Assigns OCT score to patient based on RNFLT, MRW, and GCLT
#' @param df.results formatted results to be used
#' @return none
#' @export
scoreOct = function(df.results = df.best_match)
{
  df.oct_scores = data.frame("OCT.Score" = 1:nrow(df.results))

  for(pat in 1:nrow(df.results)){
    oct.score = 0
    rnfl.count = 0
    mrw.count = 0
    gcl.count = 0

    #count score of RNFL
    rnfl.range = which(colnames(df.results) == "RNFLClass_G") : (which(colnames(df.results) == "RNFLClass_G") + 6)
    for(rnfl_val in df.results[pat,rnfl.range]){
      if(rnfl_val == "ONL"){
        rnfl.count = rnfl.count + 1
      }
    }

    #count score of MRW
    mrw.range = which(colnames(df.results) == "MRW.P.Global") : (which(colnames(df.results) == "MRW.P.Global") + 6)
    for(mrw_val in df.results[pat,mrw.range]){
      if(as.double(mrw_val) < 0.01){
        mrw.count = mrw.count + 1
      }
    }

    #count score of GCL
    gcl.range = which(colnames(df.results) == "GCL.ABN.Global") : (which(colnames(df.results) == "GCL.ABN.Global") + 6)
    for(gcl_val in df.results[pat,gcl.range]){
      if(gcl_val == T){
        gcl.count = gcl.count + 1
      }
    }

    if(rnfl.count == 1){
      oct.score = oct.score + 1
    }
    else if(rnfl.count > 2){
      oct.score = oct.score + 2
    }
    if(mrw.count == 1){
      oct.score = oct.score + 1
    }
    else if(mrw.count > 2){
      oct.score = oct.score + 2
    }
    if(gcl.count == 1){
      oct.score = oct.score + 1
    }
    else if(gcl.count > 2){
      oct.score = oct.score + 2
    }

    df.oct_scores[pat,"OCT.Score"] = oct.score
  }
  assign("df.best_match", data.frame(df.results[,1:21], df.oct_scores, df.results[,22:39]), envir = .GlobalEnv)
}

#' Reproducibility
#'
#' Assigns OCT score to patient based on RNFLT, MRW, and GCLT
#' @param df.results formatted results to be used
#' @return none
#' @export
checkCriteriaReproducibility = function()
{
  names_col = c("HAP2", "HAP2.T", "UKGTS", "UKGTS.T", "GHT", "GHT.T", "FOST", "FOST.T",
                "LOGTS", "LOGTS.T", "EAGLE", "EAGLE.T", "AGIS", "AGIS.T")
  names_row = c(">10%", "2-10%", "0.5-2%", "<0.5%")
  df.criteria_reproducibility_md = data.frame(matrix(0, nrow=length(names_row), ncol=length(names_col)))
  colnames(df.criteria_reproducibility_md) = names_col
  rownames(df.criteria_reproducibility_md) = names_row
  df.criteria_reproducibility_md = data.frame("MD"=rownames(df.criteria_reproducibility_md), df.criteria_reproducibility_md)

  names_row = c("0", "1-3", "4-5", "6")
  df.criteria_reproducibility_os = data.frame(matrix(0, nrow=length(names_row), ncol=length(names_col)))
  colnames(df.criteria_reproducibility_os) = names_col
  rownames(df.criteria_reproducibility_os) = names_row
  df.criteria_reproducibility_os = data.frame("OCT.Score"=rownames(df.criteria_reproducibility_os), df.criteria_reproducibility_os)

  followup_times = vector()

  for(pat.i in 1:nrow(df.best_match)){
    row = row2 = ""
    id = df.best_match[pat.i,"Patient.ID"]
    eye = df.best_match[pat.i,"Eye.Vf"]
    date1 = df.best_match[pat.i,"Date.Vf"]
    #print(pat.i)

    for(pat.j in which(as.integer(df.criteria_results[,"Patient.ID"]) == id)){
      date2 = as.Date.POSIXct(as.double(df.criteria_results[pat.j,"Date.Time"]))
      eye2 = df.criteria_results[pat.j,"Eye"]
      #print(pat.j)

      if((date2 > date1) && (eye2 == eye)){

        md = df.best_match[pat.i,"MD.pval"]

        if(is.na(md)){
          row = ">10%"
        }
        else if(md == 0.005){
          row = "<0.5%"
        }
        else if(md == 0.01){
          row = "0.5-2%"
        }
        else if(md == 0.02){
          row = "0.5-2%"
        }
        else if(md == 0.05){
          row = "2-10%"
        }
        else if(md == 0.1){
          row = "2-10%"
        }
        else
          print("catch")

        os = df.best_match[pat.i,"OCT.Score"]

        if(os == 0){
          row2 = "0"
        }
        else if(inside.range(os, c(1,3))){
          row2 = "1-3"
        }
        else if(inside.range(os, c(4,5))){
          row2 = "4-5"
        }
        else if(os == 6){
          row2 = "6"
        }

        if(df.best_match[pat.i,"HAP2"] == T){
          df.criteria_reproducibility_md[row,"HAP2.T"] = df.criteria_reproducibility_md[row,"HAP2.T"] + 1
          df.criteria_reproducibility_os[row2,"HAP2.T"] = df.criteria_reproducibility_os[row2,"HAP2.T"] + 1

          if(df.criteria_results[pat.j,"HAP2"] == F){
            df.criteria_reproducibility_md[row,"HAP2"] = df.criteria_reproducibility_md[row,"HAP2"] + 1
            df.criteria_reproducibility_os[row2,"HAP2"] = df.criteria_reproducibility_os[row2,"HAP2"] + 1
          }
        }

        if(df.best_match[pat.i,"UKGTS"] == T){
          df.criteria_reproducibility_md[row,"UKGTS.T"] = df.criteria_reproducibility_md[row,"UKGTS.T"] + 1
          df.criteria_reproducibility_os[row2,"UKGTS.T"] = df.criteria_reproducibility_os[row2,"UKGTS.T"] + 1

          if(df.criteria_results[pat.j,"UKGTS"] == F){
            df.criteria_reproducibility_md[row,"UKGTS"] = df.criteria_reproducibility_md[row,"UKGTS"] + 1
            df.criteria_reproducibility_os[row2,"UKGTS"] = df.criteria_reproducibility_os[row2,"UKGTS"] + 1
          }
        }

        if(df.best_match[pat.i,"GHT"] == T){
          df.criteria_reproducibility_md[row,"GHT.T"] = df.criteria_reproducibility_md[row,"GHT.T"] + 1
          df.criteria_reproducibility_os[row2,"GHT.T"] = df.criteria_reproducibility_os[row2,"GHT.T"] + 1

          if(df.criteria_results[pat.j,"GHT"] == F){
            df.criteria_reproducibility_md[row,"GHT"] = df.criteria_reproducibility_md[row,"GHT"] + 1
            df.criteria_reproducibility_os[row2,"GHT"] = df.criteria_reproducibility_os[row2,"GHT"] + 1
          }
        }

        if(df.best_match[pat.i,"FOST"] == T){
          df.criteria_reproducibility_md[row,"FOST.T"] = df.criteria_reproducibility_md[row,"FOST.T"] + 1
          df.criteria_reproducibility_os[row2,"FOST.T"] = df.criteria_reproducibility_os[row2,"FOST.T"] + 1

          if(df.criteria_results[pat.j,"FOST"] == F){
            df.criteria_reproducibility_md[row,"FOST"] = df.criteria_reproducibility_md[row,"FOST"] + 1
            df.criteria_reproducibility_os[row2,"FOST"] = df.criteria_reproducibility_os[row2,"FOST"] + 1
          }
        }

        if(df.best_match[pat.i,"LOGTS"] == T){
          df.criteria_reproducibility_md[row,"LOGTS.T"] = df.criteria_reproducibility_md[row,"LOGTS.T"] + 1
          df.criteria_reproducibility_os[row2,"LOGTS.T"] = df.criteria_reproducibility_os[row2,"LOGTS.T"] + 1

          if(df.criteria_results[pat.j,"LOGTS"] == F){
            df.criteria_reproducibility_md[row,"LOGTS"] = df.criteria_reproducibility_md[row,"LOGTS"] + 1
            df.criteria_reproducibility_os[row2,"LOGTS"] = df.criteria_reproducibility_os[row2,"LOGTS"] + 1
          }
        }

        if((df.best_match[pat.i,"HAP2"] == T) || (df.best_match[pat.i,"UKGTS"] == T) || (df.best_match[pat.i,"GHT"] == T) || (df.best_match[pat.i,"FOST"] == T) || (df.best_match[pat.i,"LOGTS"] == T)){
          #print(date2-date1)
          followup_times = c(followup_times, date2-date1)
        }
        #if(df.best_match[pat.i,"GHT"] == T){
        #  df.criteria_reproducibility_md[row,"GHT.T"] = df.criteria_reproducibility_md[row,"GHT.T"] + 1
        #  df.criteria_reproducibility_os[row2,"GHT.T"] = df.criteria_reproducibility_os[row2,"GHT.T"] + 1
        #
        #  if(df.criteria_results[pat.j,"GHT"] == F){
        #    df.criteria_reproducibility_md[row,"GHT"] = df.criteria_reproducibility_md[row,"GHT"] + 1
        #    df.criteria_reproducibility_os[row2,"GHT"] = df.criteria_reproducibility_os[row2,"GHT"] + 1
        #  }
        #}
#
        #if(df.best_match[pat.i,"GHT"] == T){
        #  df.criteria_reproducibility_md[row,"GHT.T"] = df.criteria_reproducibility_md[row,"GHT.T"] + 1
        #  df.criteria_reproducibility_os[row2,"GHT.T"] = df.criteria_reproducibility_os[row2,"GHT.T"] + 1
        #
        #  if(df.criteria_results[pat.j,"GHT"] == F){
        #    df.criteria_reproducibility_md[row,"GHT"] = df.criteria_reproducibility_md[row,"GHT"] + 1
        #    df.criteria_reproducibility_os[row2,"GHT"] = df.criteria_reproducibility_os[row2,"GHT"] + 1
        #  }
        #}
        break
      }
    }
  }

  assign("df.criteria_reproducibility_md", df.criteria_reproducibility_md, envir = .GlobalEnv)
  assign("df.criteria_reproducibility_os", df.criteria_reproducibility_os, envir = .GlobalEnv)
  assign("followup_times", followup_times, envir = .GlobalEnv)
  print(df.criteria_reproducibility_md)
  print(df.criteria_reproducibility_os)
}

age = vector()
gender = vector()
for(id in df.best_match$Patient.ID){
  for(j in 1:nrow(df.mrw_data)){
    if(id == df.mrw_data[j,"Patient.ID"]){
      age = c(age,df.mrw_data[j,"Age"])
      gender = c(gender,df.mrw_data[j,"Gender"])
      break
    }
  }
}
