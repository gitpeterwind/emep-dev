## Source code

The basic workflow is

- select tag to release
- update release/copyright notice
- update source code
- commit, tag and push

```bash
#Make sure you're on a reasonable git version (on nebula at least)
module avail | grep git
module load <git version>

## On the research repository
# select tag to release
TAG=rv4_32

# create the tag if it does not exists
#git tag $TAG cce7923 -m "EMEP course 2019"

# switch to the master branch
git checkout master

# merge new changes into master
git merge $TAG

# resolve conflicts with
git status                            # list cleared and problem files
git diff ProblemFile.f90              # inspect changes
git checkout --theirs ProblemFile.f90 # accept all changes from $TAG
git add ProblemFile.f90               # mark the file as merged
git commit                            # when no more problem files

# cherry-pick additional commits
git cherry-pick ccc19da
git cherry-pick 23038e1
git cherry-pick 3e5e922

# update release script
sed -i  "s:vers=.*\;year=....:vers=${TAG//_/.};year=`date +%Y`:" mk.OpenSource

# select chem
# e.g., EmChem16x, EmChem19a, EmChem19p...
./mk.ChemFiles EmChem19p

# create release code using four CPUs (note: don't use all CPUs - you're sharing the login node with others)
./mk.OpenSource -j4

# commit changes, tag and push
git add Makefile.Open mk.OpenSource
git commit -m "OpenSource $TAG"
git tag -a $TAG\os -m "OpenSource ${TAG//_/.} (`date +%Y%m`)"
git push --follow-tags

# replicate changes on the dev branch
git checkout dev
git cherry-pick $TAG\os
git push --follow-tags

## On the open source repository
# change to the open source repository
cd ../EMEP_MSC-W_model.${TAG//_/.}.OpenSource/code

# remove old source files
for F in *.f90; do
  grep -q ${F/.f90/.o} Makefile.SRCS || rm -i $F
done
for F in *.inc; do
  grep -q $F *.f90 || rm -i $F
done
rm -i BiomassBurningMapping_CM.inc

# make sure it compiles
make && make clean

# review and commit changes
git status
git add --all
git commit -m "OpenSource ${TAG//_/.}"
git tag ${TAG//./_} -m "OpenSource ${TAG//_/.} (`date +%Y%m`)"
git push --follow-tags
```

## Meteorology

## Other input

## Model runs