library(rsconnect)
library(BiocManager)
options(repos = BiocManager::repositories())
# Ask Adam for the token and secret
rsconnect::setAccountInfo(name='ucsdccbb',
                          token='136C7CD12693C573FB7E0E59C3DCA39B',
                          secret='yu4WMCcr5sSqAgCqiC9zEaTg+uNC96apBdbdrRp8')
options(rsconnect.max.bundle.size = 5000000000)
# getOption('rsconnect.max.bundle.size') # Check bundle size
deployApp("/Users/dchilinfuentes/CCBB_projects/shinyGWAS_Repo copy/test_prototypes//20220503_Edland_HAAS_GWAS/lewyx8-phenotype/app/", 
          appName = "20220503_Edland_HAAS_GWAS-lewyx8-phenotype", 
          account = "ucsdccbb")

