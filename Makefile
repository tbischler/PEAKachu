# Link libs for development - no needed once the lib folder is in the
# PYTHONPATH.
link_lib:
	cd bin && ln -sf ../peakachulib && cd ..
	#cd tests && ln -sf ../peakachulib && cd ..

readme_html:
	pandoc --from=markdown --to=html README.md -o README.html

readme_rst:
	grep -v "^\[!" README.md | sed -e "1d" > README.md.tmp
	pandoc --from=markdown --to=rst README.md.tmp -o README.rst
	rm README.md.tmp

readme_clean:
	rm -f README.tex README.html README.rst
	rm -f README.tex README.html README.txt

pylint:
	pylint bin/peakachu peakachulib/* tests/*

new_release:
	@echo "* Create/checkout a release branch"
	@echo "  git branch release_v0.2.X"
	@echo "  git checkout release_v0.2.X"
	@echo "* Change bin/peakachu"
	@echo "* Change setup.py"
	@echo "* Change docs/source/conf.py"
	@echo "* Change CHANGELOG.txt"
	@echo "* Create new docs"
	@echo "* Test package creation"
	@echo "* Test doc creation"
	@echo "* make package_to_pypi"
	@echo "* git add CHANGELOG.txt bin/peakachu docs/source/conf.py setup.py"
	@echo "* Commit changes e.g. 'git commit -m \"Set version to 0.2.X\"'"
	@echo "* Tag the commit e.g. 'git tag -a v0.1.X -m \"version v0.2.X\"'"
	@echo "* Merge release into dev and master"
	@echo "* Push it to github: git push"
	@echo "* Generate a new release based on this tag at"
	@echo "  https://github.com/konrad/PEAKachu/releases/new"
	@echo "* Upload new docs using 'make upload_doc'"

build:
	python3 setup.py bdist

package:
	python3 setup.py sdist
	rm -rf PEAKachu.egg-info
	ls dist/*

test:
	python3 setup.py build_ext --inplace
	python3 tests/run_all_tests.py

package_to_pypi:
	python3 setup.py sdist upload
	@echo "Go to https://pypi.python.org/pypi/PEAKachu/"

html_doc:
	cd docs && make html && cd ..

upload_doc:
	cd docs/build/html/ && zip -r PEAKachu_docs.zip * && cd ../../.. && mv docs/build/html/PEAKachu_docs.zip .
	@echo ""
	@echo "Upload PEAKachu_docs.zip at https://pypi.python.org/pypi?%3Aaction=pkg_edit&name=PEAKachu"

show_html_docs:
	firefox docs/build/html/index.html &
