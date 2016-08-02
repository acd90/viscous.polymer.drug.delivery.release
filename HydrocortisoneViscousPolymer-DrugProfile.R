# ------------------------------------------------------------------------
#                         File import
# ------------------------------------------------------------------------

#"/Users/adjuhadi/Google Drive/von Recum Lab Research/Hydrocortisone-Viscous Polymer-Concentration Release"

source(file = "./DDelivery_vfinal copy.R")
H.Raw.Data <- import_96wellplates("./", pattern = "250")

# ------------------------------------------------------------------------
#                         Annotation
# ------------------------------------------------------------------------

H.Polymers.Type <- annotate96(wellgrid = c("A1:A12, E1:E12",
                                         "B1:B12, F1:F12", 
                                         "C1:C12, G1:G12,", 
                                         "D1:D12, H1:H12"),
                            values = c("Dextran", "Alpha", "Beta", "Gamma"),
                            varname = "H.Polymers.Type")

H.Polymer.Concentrations <- annotate96(wellgrid = c("A1:H3, A7:H9", 
                                                 "A4:H6, A10:H12"),
                                    values = c("1.00mg/mL (Low)", 
                                               "100.00mg/mL (High)"),
                                    varname = "H.Polymer.Concentrations")

H.Well.Content <- annotate96(wellgrid = c("B1:D12, F1:H12", 
                                        "A1:A12, E1:E12"),
                           values = c("Experimental", "Controls"),
                           varname = "H.Well.Content")

Day.Plate <- annotate96(wellgrid = c("A1:D6", "A7:D12", "E1:H6", "E7:H12"),
                        values = c(1,2,3,4),
                        varname = "Day.Plate")

H.N <- annotate96(wellgrid = c("A1:H1, A4:H4, A7:H7, A10:H10",
                             "A2:H2, A5:H5, A8:H8, A11:H11",
                             "A3:H3, A6:H6, A9:H9, A12:H12"),
                values = c(1,2,3),
                varname = "H.N")

H.Vial.Number <- annotate96(wellgrid = c("A1:A1, A7:A7, E1:E1, E7:E7", 
                                       "A2:A2, A8:A8, E2:E2, E8:E8", 
                                       "A3:A3, A9:A9, E3:E3, E9:E9",
                                       "B1:B1, B7:B7, F1:F1, F7:F7", 
                                       "B2:B2, B8:B8, F2:F2, F8:F8",
                                       "B3:B3, B9:B9, F3:F3, F9:F9",
                                       "C1:C1, C7:C7, G1:G1, G7:G7",
                                       "C2:C2, C8:C8, G2:G2, G8:G8",
                                       "C3:C3, C9:C9, G3:G3, G9:G9",
                                       "D1:D1, D7:D7, H1:H1, H7:H7",
                                       "D2:D2, D8:D8, H2:H2, H8:H8",
                                       "D3:D3, D9:D9, H3:H3, H9:H9",
                                       "A4:A4, A10:A10, E4:E4, E10:E10",
                                       "A5:A5, A11:A11, E5:E5, E11:E11",
                                       "A6:A6, A12:A12, E6:E6, E12:E12",
                                       "B4:B4, B10:B10, F4:F4, F10:F10",
                                       "B5:B5, B11:B11, F5:F5, F11:F11",
                                       "B6:B6, B12:B12, F6:F6, F12:F12",
                                       "C4:C4, C10:C10, G4:G4, G10:G10",
                                       "C5:C5, C11:C11, G5:G5, G11:G11",
                                       "C6:C6, C12:C12, G6:G6, G12:G12",
                                       "D4:D4, D10:D10, H4:H4, H10:H10",
                                       "D5:D5, D11:D11, H5:H5, H11:H11",
                                       "D6:D6, D12:D12, H6:H6, H12:H12"),
                          values = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24),
                          varname = "H.Vial.Number")

H.Annotate96d.df <- multi_join(H.Vial.Number, 
                             H.N,
                             H.Polymers.Type, 
                             H.Polymer.Concentrations, 
                             H.Well.Content, 
                             Day.Plate, 
                             H.Raw.Data,
                             by = "well") %>% 
  get_time(., pattern = "(?<=D)[0-9]+\\.*[0-9]*") 


H.Experimental.df <- filter(H.Annotate96d.df, H.Well.Content == "Experimental")

H.Control.df <- filter(H.Annotate96d.df, Well.Content == "Controls")
# ------------------------------------------------------------------------------
#                      Calculations
# ------------------------------------------------------------------------------
#H.Annotate96d.df %>% group_by(H.Polymers.Type, Concentration,  N) %>% arrange(N)

H.Annotate96d.df <- H.Annotate96d.df %>% group_by(H.Vial.Number) %>% arrange(H.Vial.Number)

#H.Release.df <- H.Annotate96d.df %>% group_by(H.Polymers.Type, H.Polymer.Concentration, N) %>% arrange(N) %>% drug_release(., intercept = 0.03353,slope = 0.01088, vial.volume = 15, sample.volume = 15)

H.Release.df <- H.Annotate96d.df %>% group_by(H.Vial.Number) %>% arrange(Time) %>% drug_release(., intercept = 0.07148,slope = 0.02028, vial.volume = 15, sample.volume = 15)

H.Release.df$cumulative.release<-(H.Release.df$cumulative.release)/1000

Stat.Cumulative.H.Release.df<-aggregate(H.Release.df$cumulative.release, by=list(H.Release.df$Time, H.Release.df$H.Polymer.Concentrations, H.Release.df$H.Polymers.Type), FUN=mean) 
H.SD.Cumulative.Release <- aggregate(H.Release.df$cumulative.release, by=list(H.Release.df$Time, H.Release.df$H.Polymer.Concentrations, H.Release.df$H.Polymers.Type), FUN= sd )
Stat.Cumulative.H.Release.df[5] <- H.SD.Cumulative.Release[4]
colnames(Stat.Cumulative.H.Release.df) <- c("Time", "H.Polymer.Concentration", "H.Polymers.Type", "H.Mean.Cumulative.Release", "H.SD.Cumulative.Release")

#Stats.Cumulative.H.Release.df <- data.frame(row.names = Time)
#Stats <- delivery_stats(H.Release.df, H.Release.df$cumulative.release)

# ------------------------------------------------------------------------------
#                     Plotting
# ------------------------------------------------------------------------------

library(ggplot2)
require(stats)


require(ggplot2)
H.Cumulative.Release.Plot <- 
  H.Release.df$H.Polymers.Type <- factor(H.Release.df$H.Polymers.Type, levels = c("Alpha", "Beta", "Gamma", "Dextran"))
  H.Trials <- factor(H.Release.df$H.N)
  ggplot(data = H.Release.df, aes(x = Time, y=cumulative.release, shape = H.Trials))+ #, colour= H.Release.df$H.Polymer.Concentrations, shape = H.Trials))+
  geom_point(aes(color=H.Polymer.Concentrations))+
  facet_grid(.~H.Polymers.Type)+
  scale_colour_manual(values = c("deepskyblue","royalblue"))+
  #geom_smooth(se = FALSE, formula = y ~x)+
  xlab("Time (days)")+
  ylab("Cumulative Release (mg of hydrocortisone)")+
  theme_bw()+
  theme(legend.position = "bottom")+
  guides(shape = "none")
  ggsave("Cumulative.Release.Plot.Hydrocortisone.png", units = "in", width = 8, height = 5)

Stat.Cumulative.H.Release.df$H.Polymers.Type <- factor(Stat.Cumulative.H.Release.df$H.Polymers.Type, levels = c("Alpha", "Beta", "Gamma", "Dextran"))

H.Cumulative.Release.Avg.Plot <-
  ggplot(data = Stat.Cumulative.H.Release.df, aes(x = Time, y = H.Mean.Cumulative.Release, color = H.Polymer.Concentration))+
  scale_fill_brewer()+
  geom_point()+
  geom_errorbar(aes(ymin = H.Mean.Cumulative.Release-H.SD.Cumulative.Release, ymax = H.Mean.Cumulative.Release + H.SD.Cumulative.Release))+
  facet_grid(.~H.Polymers.Type)+
  scale_colour_manual(values = c("deepskyblue","royalblue"))+
  xlab("Time (days)") +
  ylab("Average Cumulative Release(mg of hydrocortisone)")+
  theme_bw()+
  theme(legend.position="bottom")+
  ggsave("Avg.Cumulative.Release.Plot.Hydrocortisone.png", width = 8, height = 5, units = "in")





