#!/bin/bash

if [ "$TRAVIS_REPO_SLUG" == "gerbaudo/SusyntHlfv" ] && [ "$TRAVIS_PULL_REQUEST" == "false" ] && [ "$TRAVIS_BRANCH" == "master" ]; then

  echo -e "Publishing doxygen...\n"
  doxygen  doc/doxygen.conf
  cp -R doc/doxygen/html $HOME/doxygen-latest

  cd $HOME
  git config --global user.email "travis@travis-ci.org"
  git config --global user.name "travis-ci"
  git clone --quiet --branch=gh-pages https://${GH_TOKEN}@github.com/gerbaudo/SusyntHlfv gh-pages > /dev/null

  cd gh-pages
  git rm -rf ./doxygen-html
  cp -Rf $HOME/doxygen-latest ./doxygen-html
  git add -f .
  git commit -m "Lastest doc on successful travis build $TRAVIS_BUILD_NUMBER auto-pushed to gh-pages"
  git push -fq origin gh-pages > /dev/null

  echo -e "Published doxygen doc to gh-pages.\n"
  
fi
