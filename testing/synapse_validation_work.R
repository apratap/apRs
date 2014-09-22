source("http://depot.sagebase.org/CRAN.R")
pkgInstall("synapseClient", "staging")

library("synapseClient")
synapseLogin()

#create a new project
p1 <- Project(name='test_apratap')
p1 <- synStore(p1)

#add annotations to a project
annotationValues(p1)
annotationValues(p1) <- list('foo'='bar')
p1 <- synStore(p1)


#test recreate the same project and add more annotations
# expected behav: the project should retain the old info
p2 <- Project(name='test_apratap')
p2 <- synStore(p2)
annotationValues(p2) <- list('test'='result')
p2 <- synStore(p2)