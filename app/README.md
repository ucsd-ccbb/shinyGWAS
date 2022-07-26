# READ ME 

This is the top directory to run the shiny app. 

```app/``` includes two subdirectories:
1. ```data/``` to store all the input files
2.  ```scripts/``` which contains all the R scripts .

and two R scripts:
1. app.R
2. shinyio_deploy.R
DO NOT remove these 2 files. These 2 files must be inside the app/ directory in order to run properly in shiny.IO


## Deploy app to shinyIO:
1. Open and edit ```app/shinyio_deploy.R``` script.
2. Change the *path* and the *appName* to in line 10 under ```deployApp()```



## Run app on RStudio:
1. Download/clone repo
2. In RStudio, open the directory ```shinyGWAS/```  in a new project session.
3. Press ```Run App``` to run the app
