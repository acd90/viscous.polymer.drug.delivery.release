# ------------------------------------------------------------------------
#                         File import
# ------------------------------------------------------------------------

# "/Users/adjuhadi/Google Drive/von Recum Lab Research/Viscous Polymer Concentrations Drug Release"

source(file = "./DDelivery_vfinal.R")
Raw.Data <- import_96wellplates("./", pattern = "235")
# ------------------------------------------------------------------------
#                         Annotation
# ------------------------------------------------------------------------

Polymers.Type <- annotate96(wellgrid = c("A1:A12, E1:E12",
                                    "B1:B12, F1:F12", 
                                    "C1:C12, G1:G12,", 
                                    "D1:D12, H1:H12"),
                       values = c("Dextran", "Alpha", "Beta", "Gamma"),
                       varname = "Polymers.Type")

Polymer.Concentration <- annotate96(wellgrid = c("A1:H3, A7:H9", 
                                            "A4:H6, A10:H12"),
                          values = c("1.00mg/mL (Low)", 
                                     "100.00mg/mL (High)"),
                          varname = "Polymer.Concentration")

Well.Content <- annotate96(wellgrid = c("B1:D12, F1:H12", 
                                       "A1:A12, E1:E12"),
                          values = c("Experimental", "Controls"),
                          varname = "Well.Content")

Day.Plate <- annotate96(wellgrid = c("A1:D6", "A7:D12", "E1:H6", "E7:H12"),
                       values = c(1,2,3,4),
                       varname = "Day.Plate")

N <- annotate96(wellgrid = c("A1:H1, A4:H4, A7:H7, A10:H10",
                             "A2:H2, A5:H5, A8:H8, A11:H11",
                             "A3:H3, A6:H6, A9:H9, A12:H12"),
                values = c(1,2,3),
                varname = "N")

Vial.Number <- annotate96(wellgrid = c("A1:A1, A7:A7, E1:E1, E7:E7", 
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
                varname = "Vial.Number")

Annotate96d.df <- multi_join(Vial.Number, 
                             N,
                             Polymers.Type, 
                             Polymer.Concentration, 
                             Well.Content, 
                             Day.Plate, 
                             Raw.Data,
                             by = "well") %>%
get_time(., pattern = "(?<=D)[0-9]+\\.*[0-9]*")

Experimental.df <- filter(Annotate96d.df, Well.Content == "Experimental")

Control.df <- filter(Annotate96d.df, Well.Content == "Controls")
# ------------------------------------------------------------------------------
#                      Calculations
# ------------------------------------------------------------------------------
# Annotate96d.df %>% group_by(Polymers.Type, Concentration,  N) %>% arrange(N)

Annotate96d.df <- Annotate96d.df %>% group_by(Vial.Number) %>% arrange(Vial.Number)

#Release2.df <- Annotate96d.df %>% group_by(Polymers.Type, Polymer.Concentration, N) %>% arrange(N) %>% drug_release(., intercept = 0.03353,slope = 0.01088, vial.volume = 15, sample.volume = 15)

Release2.df <- Annotate96d.df %>% group_by(Vial.Number) %>% arrange(Time) %>% drug_release(., intercept = 0.03353,slope = 0.01088, vial.volume = 15, sample.volume = 15)

Release2.df$cumulative.release<-(Release2.df$cumulative.release)/1000

Stat.Cumulative.Release.df<-aggregate(Release2.df$cumulative.release, by=list(Release2.df$Time, Release2.df$Polymer.Concentration, Release2.df$Polymers.Type), FUN=mean) 
SD.Cumulative.Release <- aggregate(Release2.df$cumulative.release, by=list(Release2.df$Time, Release2.df$Polymer.Concentration, Release2.df$Polymers.Type), FUN= sd )
Stat.Cumulative.Release.df[5] <- SD.Cumulative.Release[4]
colnames(Stat.Cumulative.Release.df) <- c("Time", "Polymer.Concentration", "Polymers.Type", "Mean.Cumulative.Release", "SD.Cumulative.Release")

D.Release.df<-merge(Release2.df, Stat.Cumulative.Release.df, by = c("Time", "Polymer.Concentration", "Polymers.Type"))
#Stats.Cumulative.Release.df <- data.frame(row.names = Time)
#Stats <- delivery_stats(Release2.df, Release2.df$cumulative.release)

# ------------------------------------------------------------------------------
#                     Standard Curve-drug,polymer, pbs-sds take 1
# ------------------------------------------------------------------------------

Dextran.df <- filter(Annotate96d.df, Polymers.Type == "Dextran")
Alpha.df <- filter(Annotate96d.df, Polymers.Type == "Alpha")
Beta.df <- filter(Annotate96d.df, Polymers.Type == "Beta")
Gamma.df <- filter(Annotate96d.df, Polymers.Type == "Gamma")

SCurve.Dextran.df <- Dextran.df %>% group_by(Vial.Number) %>% arrange(Time) %>% drug_release(., intercept = .16053, slope = 0.01528, vial.volume = 15, sample.volume = 15)

SCurve.Alpha.df <- Alpha.df %>% group_by(Vial.Number) %>% arrange(Time) %>% drug_release(., intercept = 0.07897,slope = 0.01361, vial.volume = 15, sample.volume = 15) 

SCurve.Beta.df <- Beta.df %>% group_by(Vial.Number) %>% arrange(Time) %>% drug_release(., intercept = 0.06339,slope = 0.01149, vial.volume = 15, sample.volume = 15) 

SCurve.Gamma.df <- Gamma.df %>% group_by(Vial.Number) %>% arrange(Time) %>% drug_release(., intercept = 0.07048,slope = 0.01116, vial.volume = 15, sample.volume = 15) 

Release3.SCurve.df <- rbind(SCurve.Dextran.df, SCurve.Alpha.df, SCurve.Beta.df, SCurve.Gamma.df)

Release3.SCurve.df$cumulative.release<-(Release3.SCurve.df$cumulative.release)/1000

SCurve.Stat.Cumulative.Release.df<-aggregate(Release3.SCurve.df$cumulative.release, by=list(Release3.SCurve.df$Time, Release3.SCurve.df$Polymer.Concentration, Release3.SCurve.df$Polymers.Type), FUN=mean) 
SCurve.SD.Cumulative.Release <- aggregate(Release3.SCurve.df$cumulative.release, by=list(Release3.SCurve.df$Time, Release3.SCurve.df$Polymer.Concentration, Release3.SCurve.df$Polymers.Type), FUN= sd )
SCurve.Stat.Cumulative.Release.df[5] <- SCurve.SD.Cumulative.Release[4]
colnames(SCurve.Stat.Cumulative.Release.df) <- c("Time", "Polymer.Concentration", "Polymers.Type", "Mean.Cumulative.Release", "SD.Cumulative.Release")

# ---------------------------------------------------------------------------
#                    Standard Curve-drug,polymer, pbs-sds take 2
# ---------------------------------------------------------------------------

SCurve2.Dextran.df <- Dextran.df %>% group_by(Vial.Number) %>% arrange(Time) %>% drug_release(., intercept = .1056, slope = 0.01499, vial.volume = 15, sample.volume = 15)

SCurve2.Alpha.df <- Alpha.df %>% group_by(Vial.Number) %>% arrange(Time) %>% drug_release(., intercept = .05906,slope = .0137, vial.volume = 15, sample.volume = 15) 

SCurve2.Beta.df <- Beta.df %>% group_by(Vial.Number) %>% arrange(Time) %>% drug_release(., intercept = .05621,slope = .01108, vial.volume = 15, sample.volume = 15) 

SCurve2.Gamma.df <- Gamma.df %>% group_by(Vial.Number) %>% arrange(Time) %>% drug_release(., intercept = .07193,slope = .0107, vial.volume = 15, sample.volume = 15) 

Release4.SCurve2.df <- rbind(SCurve2.Dextran.df, SCurve2.Alpha.df, SCurve2.Beta.df, SCurve2.Gamma.df)

Release4.SCurve2.df$cumulative.release<-(Release4.SCurve2.df$cumulative.release)/1000

SCurve2.Stat.Cumulative.Release.df<-aggregate(Release4.SCurve2.df$cumulative.release, by=list(Release4.SCurve2.df$Time, Release4.SCurve2.df$Polymer.Concentration, Release4.SCurve2.df$Polymers.Type), FUN=mean) 

SCurve2.SD.Cumulative.Release <- aggregate(Release4.SCurve2.df$cumulative.release, by=list(Release4.SCurve2.df$Time, Release4.SCurve2.df$Polymer.Concentration, Release4.SCurve2.df$Polymers.Type), FUN= sd )

SCurve2.Stat.Cumulative.Release.df[5] <- SCurve2.SD.Cumulative.Release[4]

colnames(SCurve2.Stat.Cumulative.Release.df) <- c("Time", "Polymer.Concentration", "Polymers.Type", "Mean.Cumulative.Release", "SD.Cumulative.Release")
# ------------------------------------------------------------------------------
#                     Plotting
# ------------------------------------------------------------------------------

library(ggplot2)
require(stats)

#quartzFont(helvectica.neue = "helvectica neue")

#Release2.df$Polymer.Type <- factor(Release2.df$Polymers.Type, ordered=as.ordered("Dextran", "Alpha","Beta","Gamma"))

Release2.df$Polymers.Type <- factor(Release2.df$Polymers.Type, levels = c("Alpha", "Beta", "Gamma", "Dextran"))
Release2.df$Drug = "Dexamethasone"
require(ggplot2)
Cumulative.Release.Plot <- 
  Trials <- factor(Release2.df$N)
  ggplot(data = Release2.df, aes(x = Time, y = cumulative.release, shape = Trials, color = Polymer.Concentration))+ 
  geom_point()+
  facet_grid(.~Polymers.Type)+
  #scale_colour_manual(values =c("deepskyblue","royalblue"))+
  #facet_wrap()
  #geom_smooth(se = FALSE, formula = y ~x)+
  xlab("Time (days)")+
  ylab("Cumulative Release (mg of dexamethasone)")+
  theme_bw()+
  theme(legend.position = "bottom")+
  guides(shape = "none")
  ggsave("Cumulative.Release.Plot.Dexamethasone.png", units = "in", width = 8, height = 5)

  Trials <- factor(Release2.df$N)
    geom_point()+
#___________________________________________________________
#Create 4 different graphs into 1
Release.Dextran.df <- filter(Release2.df, Polymers.Type == "Dextran")
Release.Alpha.df <- filter(Release2.df, Polymers.Type == "Alpha")
Release.Beta.df <- filter(Release2.df, Polymers.Type == "Beta")
Release.Gamma.df <- filter(Release2.df, Polymers.Type == "Gamma")

#Dextran  
Release.Dextran <- 
  D.Trials <- factor(Release.Dextran.df$N)
Release.Dextran <-
  ggplot(data = Release.Dextran.df, aes(x = Time, y = cumulative.release, color =Polymer.Concentration, shape= D.Trials))+
  geom_point()+
  scale_colour_manual(values = c("red","red4"))+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "vertical", axis.title.y=element_blank(), legend.key.height = unit(0.05,"in"))+
  coord_cartesian(ylim=c(0,4.1))+
  xlab("Time (days)")+
  guides(shape = "none")+
  ggtitle("Dextran")

#Alpha
  Release.Alpha <- 
    A.Trials <- factor(Release.Alpha.df$N)
  Release.Alpha <-
  ggplot(data = Release.Alpha.df, aes(x = Time, y = cumulative.release, color =Polymer.Concentration, shape= A.Trials))+
    geom_point()+
    scale_colour_manual(values = c("purple","darkorchid4"))+
    theme_bw()+
    xlab("Time (days)")+
    coord_cartesian(ylim = c(0,4.1))+
    ylab("Cumulative Release (mg of dexamethasone)")+
    ggtitle("Alpha")+
    theme(legend.position = "bottom", legend.direction = "vertical", axis.title.y=element_blank(), legend.key.height = unit(0.05,"in"), axis.title.y	= )+
    guides(shape = "none")
   

  #Beta
  Release.Beta <- 
    B.Trials <- factor(Release.Beta.df$N)
  Release.Beta <-
  ggplot(data = Release.Beta.df, aes(x = Time, y = cumulative.release, color =Polymer.Concentration, shape= B.Trials))+
    geom_point()+
    scale_colour_manual(values = c("deepskyblue","blue"))+
    theme_bw()+
    theme(legend.position = "bottom", legend.direction = "vertical", axis.title.y=element_blank(), legend.key.height = unit(0.05,"in"))+
    coord_cartesian(ylim=c(0,4.1))+
    guides(shape = "none")+
    xlab("Time (days)")+
    ggtitle("Beta")
  
#Gamma
  Release.Gamma <- 
    G.Trials <- factor(Release.Gamma.df$N)
  Release.Gamma <- ggplot(data = Release.Gamma.df, aes(x = Time, y = cumulative.release, color =Polymer.Concentration, shape= G.Trials))+
    geom_point()+
    scale_colour_manual(values = c("green2","green4"))+
    coord_cartesian(ylim=c(0,4.1))+
    theme_bw()+
    theme(legend.position = "bottom", legend.direction = "vertical", axis.title.y=element_blank(), legend.key.height = unit(0.05,"in"))+
    guides(shape = "none")+
    xlab("Time (days)")+
    ggtitle("Gamma")


#--------------------------------------------------------------------- 
  # Multiple plot function
  # ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
  # - cols:   Number of columns in layout
  # - layout: A matrix specifying the layout. If present, 'cols' is ignored.
  #
  # If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
  # then plot 1 will go in the upper left, 2 will go in the upper right, and
  # 3 will go all the way across the bottom.
  #
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  #--------------------------------------------------------------------
require(ggplot2)
require(grid)
require(gridExtra)
library(devtools)
library(easyGgplot2)
ggplot2.multiplot(Release.Alpha,Release.Beta,Release.Gamma,Release.Dextran, cols= 4)+

require(cowplot)
plot_grid(Release.Alpha, Release.Beta, Release.Gamma, Release.Dextran, ncol = 2, nrow=1)

#----------------------------------------------------------------------
  #Using different standard curves by polymer
Scurve.Cumulative.Release.Plot <- 
    Trials <- factor(Release3.SCurve.df$N)
  ggplot(data = Release3.SCurve.df, aes(x = Time, y = cumulative.release, color = Polymer.Concentration, shape = Trials))+ 
    geom_point()+
    facet_grid(.~Polymers.Type)+
    #geom_smooth(se = FALSE, formula = y ~x)+
    xlab("time (days)")+
    ylab("cumulative.release (mg)")
  
#Using different standard curves by polymer take 2
  Scurve2.Cumulative.Release.Plot <- 
    Trials <- factor(Release4.SCurve2.df$N)
  ggplot(data = Release4.SCurve2.df, aes(x = Time, y = cumulative.release, color = Polymer.Concentration, shape = Trials))+ 
    geom_point()+
    facet_grid(.~Polymers.Type)+
    #geom_smooth(se = FALSE, formula = y ~x)+
    xlab("time (days)")+
    ylab("cumulative.release (mg)")

  Stat.Cumulative.Release.df$Polymers.Type <- factor(Stat.Cumulative.Release.df$Polymers.Type, levels = c("Alpha", "Beta", "Gamma", "Dextran"))
  
Cumulative.Release.Avg.Plot <-
  ggplot(data = Stat.Cumulative.Release.df, aes(x = Time, y = Mean.Cumulative.Release, color = Polymer.Concentration))+
  scale_fill_brewer()+
  geom_point()+
  geom_errorbar(aes(ymin = Mean.Cumulative.Release-SD.Cumulative.Release, ymax = Mean.Cumulative.Release + SD.Cumulative.Release))+
  facet_grid(.~Polymers.Type)+
  scale_colour_manual(values = c("deepskyblue","royalblue"))+
  xlab("Time (days)") +
  ylab("Average Cumulative Release(mg of dexamethasone)")+
  theme_bw()+
  theme(legend.position = "bottom")
  ggsave("Avg.Cumulative.Release.Plot.Dexamethasone.png", width = 8, height = 5, units = "in")
  
  
Time.Point.Plot<- ggplot(data = Release2.df, aes(x = time, y = timepoint.release))+ 
  geom_point(aes(color = Polymer.Concentration)) +
  facet_grid(.~Polymers.Type)+
  scale_y_log10()


#, labeller = label_parsed(labels, multi_line = TRUE)) 
  #geom_smooth(mapping = NULL, data = NULL, stat = "smooth",  method = "lm", formula = y ~ x, se = TRUE, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) +
  #theme_bw()

absorbance.exp <-
  ggplot(data = experimental.df, aes(x = time, y = absorbance)) +
  geom_point(aes(color = Polymer.Concentration)) +
  facet_grid(.~Polymers.Type, labeller = label_parsed) +
  geom_smooth(mapping = NULL, data = NULL, stat = "smooth",  method = "lm", formula = y ~ x, se = TRUE, na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) +
  theme_bw()

absorbance.exp +
  geom_point(aes(color = Polymers.Type))+
  facet_grid(.~Polymer.Concentration, labeller =label_parsed) +
  theme_bw()

absorbance.ctrl <- control.df %>%
  ggplot(., aes(x = time, y = absorbance)) + 
  geom_point()

absorbance.ctrl +
  theme_bw()

#----------------------------------------------------------------------------------