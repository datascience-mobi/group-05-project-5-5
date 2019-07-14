Are differentially methylated regions within genes associated with mantle cell lymphoma?
========================================================
author: Pascal Lafrenz, Mari Hambardzumyan, Lea Herzel, Franziska Lam
date: July 24. 2019
autosize: true

Project milestones
========================================================

$$~$$




1. Data processing 
2. Data normalization and visualization 
3. Data reduction
4. Regression and interpretation

1. Data processing
========================================================
Initial goal: retail 90% of the information after data processing.
$$~$$

How many genes across 10 samples in total?











































































```
processing file: final Presentation.Rpres
Registered S3 methods overwritten by 'ggplot2':
  method         from 
  [.quosures     rlang
  c.quosures     rlang
  print.quosures rlang
-- Attaching packages --------------------------------------------------------------------------------- tidyverse 1.2.1 --
v ggplot2 3.1.1     v purrr   0.3.2
v tibble  2.1.3     v dplyr   0.8.1
v tidyr   0.8.3     v stringr 1.4.0
v readr   1.3.1     v forcats 0.4.0
-- Conflicts ------------------------------------------------------------------------------------ tidyverse_conflicts() --
x dplyr::filter()     masks stats::filter()
x dplyr::group_rows() masks kableExtra::group_rows()
x dplyr::lag()        masks stats::lag()

Attaching package: 'gridExtra'

The following object is masked from 'package:dplyr':

    combine


Attaching package: 'plotly'

The following object is masked from 'package:ggplot2':

    last_plot

The following object is masked from 'package:stats':

    filter

The following object is masked from 'package:graphics':

    layout


Attaching package: 'gplots'

The following object is masked from 'package:stats':

    lowess

corrplot 0.84 loaded
Quitting from lines 49-53 (final Presentation.Rpres) 
Fehler in gzfile(file, "rb") : kann Verbindung nicht öffnen
Ruft auf: knit ... withCallingHandlers -> withVisible -> eval -> eval -> readRDS -> gzfile
Zusätzlich: Warnmeldungen:
1: package 'tidyverse' was built under R version 3.6.1 
2: package 'ggrepel' was built under R version 3.6.1 
3: package 'gplots' was built under R version 3.6.1 
4: package 'sandwich' was built under R version 3.6.1 
5: package 'corrplot' was built under R version 3.6.1 
6: package 'jtools' was built under R version 3.6.1 
Ausführung angehalten
```
