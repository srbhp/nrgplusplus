#!/bin/bash

rm -rf "$1"

git clone https://github.com/brunowu/doxygen_sphinx_example.git $1

cd "$1" 

git checkout gh-pages

cp -rf "$2/" .

git add --all .

git commit -m "message[$((1 + $RANDOM % 100000))]: update the webpages"

git push origin gh-pages

cd ..

rm -rf "$1"
