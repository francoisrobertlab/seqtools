if grep -Fq "ROBERTF_MODULES_DIR" ~/.bash_profile ; then
  echo "Modules directory present in .bash_profile"
else
  echo "Adding modules directory to .bash_profile"
  echo '## Robert Lab Modules ##' >> ~/.bash_profile
  echo 'ROBERTF_MODULES_DIR=~/projects/def-robertf/modules' >> ~/.bash_profile
  echo 'if [ -d "$ROBERTF_MODULES_DIR" ]; then' >> ~/.bash_profile
  echo '  module use $ROBERTF_MODULES_DIR' >> ~/.bash_profile
  echo 'fi' >> ~/.bash_profile
  echo '' >> ~/.bash_profile
fi

