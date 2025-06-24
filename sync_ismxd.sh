#!/bin/bash

rsync -avz -e ssh --exclude work --exclude .git ./ ismxd:/home/25cr05tm/ActFEM/