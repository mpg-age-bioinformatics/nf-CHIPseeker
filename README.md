# nf-CHIPseeker

The current function is to annoate the peak for the output files from diffbind[nf-diffbind](https://github.com/mpg-age-bioinformatics/nf-diffbind?tab=readme-ov-file#for-annotating-the-peaks)

## Contributing

Make a commit, check the last tag, add a new one, push it and make a release:
```
git add -A . && git commit -m "<message>" && git push
git describe --abbrev=0 --tags
git tag -e -a <tag> HEAD
git push origin --tags
gh release create <tag> 
```