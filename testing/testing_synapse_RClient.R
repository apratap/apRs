#testing synapse connection

remove.packages("synapseClient")
source("http://depot.sagebase.org/CRAN.R")
pkgInstall("synapseClient",stack="staging")

sessionInfo()

library(synapseClient)


path.expand('~/.synapseConfig')
file.exists('~/.synapseConfig')
synapseClient:::ConfigParser()

library(synapseClient)
synapseLogin()

myproject = Project(name='test_apratap1')
myproject = synStore(myproject)
pid <- propertyValue(myproject,"id")








wikiPage <- WikiPage(owner=myproject,
                     title = "test_wiki",
                     markdown = "*TEST*",
                     fileHandles = list(fh$id))
wikiPage <- synStore(wikiPage)
wikiPage2 <- synGetWiki(myproject)

wikiPage2


file <- synGet('syn2176182')
evaluation<-Evaluation(name="test evaluation 7", status="PLANNED", contentSource=pid)
evaluation<-synStore(evaluation)
eid<-propertyValue(evaluation,"id")

#open the elaluation
propertyValue(evaluation, "status")<-"OPEN"
evaluation<-synStore(evaluation)

#add user to the evaluation
synapseClient:::.allowParticipation(eid,"PUBLIC")
synRestPOST(sprintf("/evaluation/%s/participant",eid),list())

submit(evaluation,file)

