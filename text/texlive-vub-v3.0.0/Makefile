ZIP=beamercolorthemevub.sty \
    beamerfontthemevub.sty \
    beamerinnerthemevub.sty \
    beamerouterthemevub.sty \
    beamerthemevub.sty \
    dept/vub-ai.pdf \
    dept/vub-brubotics-cmyk.pdf \
    dept/vub-etro-cmyk.pdf \
    dept/vub-soft.pdf \
    dept/vub-indi-cmyk.pdf \
    \
	ulb_logo_3lp.eps \
	bruface.png \
	vub_logo_cmyk.pdf \
	\
	vub.sty \
	vubprivate.sty \
	bruface.sty \
	\
	README.md \

all: texlive-vub.zip templates

texlive-vub.zip: $(ZIP)
	echo making $@
	rm -f $@
	zip $@ $(ZIP)

.PHONY: templates publish

templates:
	$(MAKE) -C templates

publish:
	$(MAKE) -C templates publish
