#!/bin/bash
set -e
rm -rf m4rie
ln -sf src m4rie
autoreconf --install
