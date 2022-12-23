TEMPLATES=vub-article \
		  vub-presentation \
		  vub-thesis \
		  bruface-thesis \

.PHONY: ../texlive-vub.zip publish

all: $(addprefix template-, $(addsuffix .zip,$(TEMPLATES)))

publish: $(addprefix .publish-, $(TEMPLATES))

../texlive-vub.zip:
	$(MAKE) -C .. texlive-vub.zip

template-%.zip: ../texlive-vub.zip %.tex
	echo making $* from $<
	cp $< $@
	mkdir -p $*
	cp $*.tex $*/main.tex
	zip $@ -j $*/main.tex
	rm $*/main.tex
	rmdir $*

.publish-%: template-%.zip
	../.ci/push-to-overleaf.sh $*
	touch $@
