# This is R code for the clustering analysis of Uvira cholera data
# 2016-2020. This script produces Figure 1 (epidemic curve).
# I combine the epicurve with the seasonality assessment graph
# Author: R Ratnayake, 2023

# Epicurve and seasonality, Uvira

# This script produces the epicurve for 2016-2020, by testing status

# Run uvira data prep script first
# Filter to 2016-2020 and replace NA with Untested
uvira.test <- uvira %>% 
  filter(year>=2016 & year<=2020)

uvira.test$result <- uvira.test$result %>% 
  replace_na('Untested')
  
result_epi_month <- incidence2::incidence(
  x = uvira.test, 
  date_index = date_adm,
  interval = "month",
  groups = result,
  na_as_group = TRUE)

summary(result_epi_month)

result_epi_month$result = factor(result_epi_month$result,
                                 levels=c("Untested", 
                                          "Negative",
                                          "Positive"))
#produce the epicurve

epicurve <- plot(result_epi_month, 
     fill = result,
     #show_cases=TRUE,
     color="black",
     alpha = 0.9,
     size = 11,
     legend = "bottom",
     date_format="%Y-%b\n(Week %W)",
     centre_dates=FALSE) +
     xlab("\nMonth of case presentation") +
     scale_y_continuous(breaks=seq(0,300,50),
                        limits=c(0,300),
                        expand=c(0,0)) +
     scale_fill_manual(values=c("#899DA4","#DC863B","#C93312"), 
                       name = "Test status") +
     theme(axis.line = element_line(colour = "black"),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           panel.background = element_blank(),
           legend.position = c(0.9, 0.9),
           axis.text = element_text(size=11))

# Combine graphs

#run seasonality analysis
tiff("epicurve_seasonality.tiff", units="in", width=11.69, height=8.27, res=300)
ggdraw() + draw_plot(epicurve) +
  draw_plot(p3, x = 0.55, y = .74, width = .25, height = .25)
dev.off()

# specify the date of the first case: "Friday 01 Jan, 2016" and breaks
format(min(uvira.test$date_adm, na.rm=T), "%A %d %b, %Y")

(monthly_breaks_centered <- seq.Date(
  from = floor_date(min(uvira.test$date_adm, na.rm=T), 
                    "week", week_start = 1), 
  to   = ceiling_date(max(uvira.test$date_adm, na.rm=T), 
                      "week", week_start = 1),
  by   = "month"))

(monthly_breaks <- seq.Date(from = as.Date("2016-01-01"),
                           to = as.Date("2020-12-31"),
                           by = "month"))

(weekly_breaks <- seq.Date(from = as.Date("2016-01-01"),
                            to = as.Date("2020-12-31"),
                            by = "week"))

# set order of test result status
uvira.test$result = factor(uvira.test$result,
                           levels=c("Untested", "Negative", "Positive"))

(epicurve2 <- 
  ggplot(uvira.test) +
  geom_histogram(aes(x=date_adm, group = result, fill = result),
                 breaks = monthly_breaks, 
                 closed = "left",
                 color = "black", size=0.1) + 
  scale_x_date(limits = as.Date(c("2016-01-01", "2020-12-31")),
               expand = c(0,0),
               date_breaks = "6 months",
               date_labels = "\n%Y\n%b") +
  scale_y_continuous(breaks=seq(0,300,50),
                     limits=c(0,300),
                     expand=c(0,0)) +
  scale_fill_manual(values=c("#899DA4","#DC863B","#C93312"), 
                    name = "Test status") +
  labs(x = "\nMonth of case presentation", y = "\nMonthly suspected cases") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.85, 0.9),
        axis.text = element_text(size=11))) 

# Combine graphs (ggplot version)
tiff("epicurve_seasonality.tiff", units="in", width=11.69, height=8.27, res=300)
ggdraw() + draw_plot(epicurve2) +
  draw_plot(p3, x = 0.50, y = .74, width = .25, height = .25)
dev.off()