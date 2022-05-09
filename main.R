library(tercenApi)
library(tercen)
library(dplyr)
library(progressr)
library("future.apply")

# options("tercen.serviceUri"="http://172.28.0.1:5400/api/v1/")
# # # http://127.0.0.1:5400/test-team/w/073510448c675ef923a0b55ca20ba1c0/ds/9fb0dd32-20d1-4daa-8701-e5766bfb425c
# options("tercen.workflowId"= "073510448c675ef923a0b55ca20ba1c0")
# options("tercen.stepId"= "9fb0dd32-20d1-4daa-8701-e5766bfb425c")

ctx = tercenCtx()

nCpus = availableCores() 
nCpusRequested = 4
ctx$requestResources(nCpus=nCpusRequested, ram=500000000, ram_per_cpu=500000000)
nCpusReceived = availableCores() 

msg = paste0("nCpus=" , nCpus , " nCpusRequested=", nCpusRequested, " nCpusReceived=", nCpusReceived)

ctx$log(msg)

ctx  %>%
  select(.y, .ci, .ri) %>%
  group_by(.ci, .ri) %>%
  summarise(mean = mean(.y)) %>%
  ctx$addNamespace() %>%
  ctx$save()
