EMAIL=$1
if [[ ! "$EMAIL" =~ ^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,4}$ ]]
then
    echo "You must supply your email address as the first argument"
    exit 1
fi

if grep -Fq "USER_EMAIL" ~/.bash_profile ; then
  echo "Email address environment variables present in .bash_profile"
else
  echo "Adding email address to environment variables in .bash_profile"
  echo "USER_EMAIL=$EMAIL" >> ~/.bash_profile
  echo "export USER_EMAIL" >> ~/.bash_profile
  echo "" >> ~/.bash_profile
fi

if grep -Fq "ROBERTF_MODULES_DIR" ~/.bash_profile ; then
  echo "Modules directory present in .bash_profile"
else
  echo "Adding modules directory to .bash_profile"
  echo "## Robert Lab Modules ##" >> ~/.bash_profile
  echo "ROBERTF_MODULES_DIR=~/projects/def-robertf/modules" >> ~/.bash_profile
  echo 'if [ -d "$ROBERTF_MODULES_DIR" ]; then' >> ~/.bash_profile
  echo '  module use $ROBERTF_MODULES_DIR' >> ~/.bash_profile
  echo 'fi' >> ~/.bash_profile
  echo "" >> ~/.bash_profile
fi

