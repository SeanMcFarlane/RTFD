#!/bin/bash

# This script deletes the project from the test environment and reuploads it.
#ssh -i ~/.ssh/id_rsa seanmcf@13.88.243.86
#rm -rf ./RTFD_2
#exit
scp -r ../RTFD_2 seanmcf@13.88.243.86:~/