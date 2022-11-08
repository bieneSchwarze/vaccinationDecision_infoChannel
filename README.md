# Impact on Information Channels on the CoV-2 Vaccination Decision (part of the Study SOEP-RKI-2)
## Data Preparation and Analysis 
for a joint project between SZinn and Susanne Jordan & Sarah Jane Böttger (RKI)

For our analysis we use data from the German RKI-SOEP-2 Project ``Corona-Monitoring bundesweit" (CoMoBu-Studie)'', see https://www.rki.de/DE/Content/Gesundheitsmonitoring/Studien/lid/lid_node.html
This data is sensitive and access has to be asked for from the RKI or the SOEP. The study design is described in: Bartig, S., Brücker, H., Butschalowsky, H., Danne, C., Gößwald, A., Goßner, L., Grabka, M., . . . , Zinn, S. (2022) Corona Monitoring Nationwide (RKI-SOEP-2): Seroepidemiological Study on the Spread of SARS-CoV-2 Across Germany. Jahrbücher für Nationalökonomie und Statistik. https://doi.org/10.1515/jbnst-2022-0047.

In this Git repository you find the R source code (R version 4.1.12) used for data preparation and model estimation. 
The main file is *information_mainfile.R*. From within this file a function *fitModel()* for estimating a survey-weighted logit regression model is called. 
The respective source code is included in the file *models.R*.
Before conducting model estimation, we impute missing values by mice (package version 3.14.0) and CART. 
Estimation results are pooled by Rubin's combining rules. 
For fitting the survey-weighted logit regression models, we use Lumley's survey package (version 4.1.1) in R.


