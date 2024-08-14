# Data: Three transects with 2 sampling depths each = six different microplastic samples
# for each sample: remove plastic pieces larger than 5 mm (not microplastics)

library(dplyr)

# Transect 1 - 0-5 cm depth -----
dat1.1.a = read.csv("RolfEtAl2022_data/imagelab_data/0-5/ML1_0-5cm_ATR1.csv", header = TRUE, stringsAsFactors = TRUE,
                    sep = ";", fileEncoding="latin1")
dat1.1.b = read.csv("RolfEtAl2022_data/imagelab_data/0-5/ML1_0-5cm_ATR2.csv", header = TRUE, stringsAsFactors = TRUE,
                   sep = ";", fileEncoding="latin1")
dat1.1 = rbind(dat1.1.a,dat1.1.b)
# remove larger plastic pieces
dat1.1 = dat1.1[dat1.1$Länge <= 5000,]
# remove particles without measurements
dat1.1 = dat1.1[!is.na(dat1.1$Form.Faser..Fragment..Kugel.),]

# get frequencies of shapes
dat1.1 %>%
  group_by(Form.Faser..Fragment..Kugel.) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# get frequencies of polymer types
dat1.1 %>%
  group_by(Analyse) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# get average and median lengths
summary(dat1.1$Länge)

# get average and median widths
summary(dat1.1$Breite)


# Transect 2 - 0-5 cm depth -----
dat2.1.a = read.csv("RolfEtAl2022_data/imagelab_data/0-5/ML2_0-5cm_ATR1.csv", header = TRUE, stringsAsFactors = TRUE,
                    sep = ";", fileEncoding="latin1")
dat2.1.b = read.csv("RolfEtAl2022_data/imagelab_data/0-5/ML2_0-5cm_ATR2.csv", header = TRUE, stringsAsFactors = TRUE,
                    sep = ";", fileEncoding="latin1")
dat2.1.c = read.csv("RolfEtAl2022_data/imagelab_data/0-5/ML2_0-5cm_ATR3.csv", header = TRUE, stringsAsFactors = TRUE,
                    sep = ";", fileEncoding="latin1")
dat2.1.d = read.csv("RolfEtAl2022_data/imagelab_data/0-5/ML2_0-5cm_ATR4.csv", header = TRUE, stringsAsFactors = TRUE,
                    sep = ";", fileEncoding="latin1")

dat2.1 = rbind(dat2.1.a,dat2.1.b, dat2.1.c, dat2.1.d)
# remove larger plastic pieces
dat2.1 = dat2.1[dat2.1$Länge <= 5000,]
# remove particles without measurements
dat2.1 = dat2.1[!is.na(dat2.1$Form.Faser..Fragment..Kugel.),]

# get frequencies of shapes
dat2.1 %>%
  group_by(Form.Faser..Fragment..Kugel.) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# get frequencies of polymer types
dat2.1 %>%
  group_by(Analyse) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# get average and median lengths
summary(dat2.1$Länge)

# get average and median widths
summary(dat2.1$Breite)


# Transect 3 - 0-5 cm depth -----
dat3.1.a = read.csv("RolfEtAl2022_data/imagelab_data/0-5/ML3_0-5cm_ATR1.csv", header = TRUE, stringsAsFactors = TRUE,
                    sep = ";", fileEncoding="latin1")
dat3.1.b = read.csv("RolfEtAl2022_data/imagelab_data/0-5/ML3_0-5cm_ATR2.csv", header = TRUE, stringsAsFactors = TRUE,
                    sep = ";", fileEncoding="latin1")

dat3.1 = rbind(dat3.1.a,dat3.1.b)
# remove larger plastic pieces
dat3.1 = dat3.1[dat3.1$Länge <= 5000,]
# remove particles without measurements
dat3.1 = dat3.1[!is.na(dat3.1$Form.Faser..Fragment..Kugel.),]

# get frequencies of shapes
dat3.1 %>%
  group_by(Form.Faser..Fragment..Kugel.) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# get frequencies of polymer types
dat3.1 %>%
  group_by(Analyse) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# get average and median lengths
summary(dat3.1$Länge)

# get average and median widths
summary(dat3.1$Breite)

# Transect 1 - 5-20 cm depth -----
dat1.2.a = read.csv("RolfEtAl2022_data/imagelab_data/5-20/ML1_5-20cm_ATR1.csv", header = TRUE, stringsAsFactors = TRUE,
                    sep = ";", fileEncoding="latin1")
dat1.2.b = read.csv("RolfEtAl2022_data/imagelab_data/5-20/ML1_5-20cm_ATR2.csv", header = TRUE, stringsAsFactors = TRUE,
                    sep = ";", fileEncoding="latin1")
dat1.2.c = read.csv("RolfEtAl2022_data/imagelab_data/5-20/ML1_5-20cm_ATR3.csv", header = TRUE, stringsAsFactors = TRUE,
                    sep = ";", fileEncoding="latin1")


dat1.2 = rbind(dat1.2.a,dat1.2.b, dat1.2.c)
# remove larger plastic pieces
dat1.2 = dat1.2[dat1.2$Länge <= 5000,]
# remove particles without measurements
dat1.2 = dat1.2[!is.na(dat1.2$Form.Faser..Fragment..Kugel.),]

# get frequencies of shapes
dat1.2 %>%
  group_by(Form.Faser..Fragment..Kugel.) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# get frequencies of polymer types
dat1.2 %>%
  group_by(Analyse) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# get average and median lengths
summary(dat1.2$Länge)

# get average and median widths
summary(dat1.2$Breite)

# Transect 2 - 5-20 cm depth -----
dat2.2.a = read.csv("RolfEtAl2022_data/imagelab_data/5-20/ML2_5-20cm_ATR1.csv", header = TRUE, stringsAsFactors = TRUE,
                    sep = ";", fileEncoding="latin1")
dat2.2.b = read.csv("RolfEtAl2022_data/imagelab_data/5-20/ML2_5-20cm_ATR2.csv", header = TRUE, stringsAsFactors = TRUE,
                    sep = ";", fileEncoding="latin1")
  
dat2.2 = rbind(dat2.2.a, dat2.2.b)

# remove larger plastic pieces
dat2.2 = dat2.2[dat2.2$Länge <= 5000,]
# remove particles without measurements
dat2.2 = dat2.2[!is.na(dat2.2$Form.Faser..Fragment..Kugel.),]

# get frequencies of shapes
dat2.2 %>%
  group_by(Form.Faser..Fragment..Kugel.) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# get frequencies of polymer types
dat2.2 %>%
  group_by(Analyse) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# get average and median lengths
summary(dat2.2$Länge)

# get average and median widths
summary(dat2.2$Breite)

# Transect 3 - 5-20 cm depth -----
dat3.2.a = read.csv("RolfEtAl2022_data/imagelab_data/5-20/ML3_5-20cm_ATR1.csv", header = TRUE, stringsAsFactors = TRUE,
                    sep = ";", fileEncoding="latin1")
dat3.2.b = read.csv("RolfEtAl2022_data/imagelab_data/5-20/ML3_5-20cm_ATR2.csv", header = TRUE, stringsAsFactors = TRUE,
                    sep = ";", fileEncoding="latin1")

dat3.2 = rbind(dat3.2.a, dat3.2.b)

# remove larger plastic pieces
dat3.2 = dat3.2[dat3.2$Länge <= 5000,]
# remove particles without measurements
dat3.2 = dat3.2[!is.na(dat3.2$Form.Faser..Fragment..Kugel.),]

# get frequencies of shapes
dat3.2 %>%
  group_by(Form.Faser..Fragment..Kugel.) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# get frequencies of polymer types
dat3.2 %>%
  group_by(Analyse) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

# get average and median lengths
summary(dat3.2$Länge)

# get average and median widths
summary(dat3.2$Breite)
