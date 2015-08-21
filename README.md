

Quickstart
----------

```sh
git clone git@github.com:susynt/susynt-read.git
cd susynt-read
gco -b n0212-cutflow n0212
git clone git@github.com:gerbaudo/Superflow.git
git clone git@github.com:gerbaudo/susynt-ss3l.git
setupATLAS
source bash/setup_area.sh
source bash/setup_release.sh

# numpy is needed for some of the scripts
localSetupSFT --cmtConfig=x86_64-slc6-gcc48-opt pyanalysis/1.4_python2.7

```

[doxygen doc](http://gerbaudo.github.io/susynt-ss3l/doxygen-html/)

davide.gerbaudo@gmail.com
Aug 2015
