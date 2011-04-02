#!/bin/tcsh

# download CLI tools
# http://aws.amazon.com/developertools/AWS-Identity-and-Access-Management/4143

setenv JAVA_HOME /usr/
setenv AWS_IAM_HOME ~/Downloads/IAMCli-1.1.0
setenv PATH $AWS_IAM_HOME/bin:$PATH
setenv AWS_CREDENTIAL_FILE /Users/depristo/Desktop/broadLocal/GATK/trunk/account-key

setenv CREATE_GROUPS false
setenv CREATE_IGV_USER false
setenv UPDATE_USER_KEYS false
setenv UPDATE_USER_POLICY true

# Create the administrators group:
# we aren't actually using this, in fact
if ( $CREATE_GROUPS == true ) then
iam-groupcreate -g Admins
iam-grouplistbypath
iam-groupuploadpolicy -g Admins -p AdminsGroupPolicy -f GroupPolicy.txt
iam-grouplistpolicies -g Admins
endif

# Create the IGV user -- uncomment if the IGV user needs to be created from scratch
# update the secret key
if $CREATE_IGV_USER == true then
iam-usercreate -u IGV -k -v > IGV_cred.txt
endif

# the user access and secret keys are in the IGV source file IGVRunReport.java
# and must be updated to be the most current ones
if $UPDATE_USER_KEYS == true then
iam-userdelkey -u IGV -k $1 # $1 -> current access key
iam-useraddkey -u IGV > IGV_cred.txt 
cat IGV_cred.txt
endif

echo "IGV user policies"
if $UPDATE_USER_POLICY == true then
echo "Deleting policy"
iam-userdelpolicy -u IGV -p IGVRunReportUploading
iam-useruploadpolicy -u IGV -p IGVRunReportUploading -f IGVPolicy.txt
endif
iam-userlistpolicies -u IGV -v 
