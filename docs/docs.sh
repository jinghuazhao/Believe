#!/usr/bin/bash

function github_pages()
{
  module load python/3.8
  source ~/rds/public_databases/software/py38/bin/activate
# module load python/3.7
# source ~/COVID-19/py37/bin/activate
# pip install mkdocs-mermaid2-plugin
}

module load ceuadmin/libssh/0.10.6-icelake
module load ceuadmin/openssh/9.7p1-icelake

github_pages
mkdocs build
mkdocs gh-deploy

git add .gitignore
git add docs
git add README.md
git add mkdocs.yml
git commit -m "Site files"
git add 0_uitls.sh 1_desc.sh 2_gcta.sb
git commit -m "script"
git push
