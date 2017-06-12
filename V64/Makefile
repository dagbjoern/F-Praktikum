ifeq (,$(shell sh -c 'cygpath --version 2> /dev/null'))
  # Unix
  pwd := $$(pwd)
  translate = $1
else
  # Windows mit MSys2/Cygwin
  pwd := $$(cygpath -m "$$(pwd)")
  translate = $(shell echo '$1' | sed 's/:/;/g')
endif

all: build/main.pdf

# hier Python-Skripte:
build/plot1.pdf: auswertung.py messungglas.txt Messwertewinkel.txt matplotlibrc header-matplotlib.tex | build
	TEXINPUTS="$(call translate,$(pwd):)" python auswertung.py

#build/plot2.pdf: build/plot1.pdf | build


#hier weitere Abhängigkeiten für build/main.pdf deklarieren:
build/main.pdf: build/plot1.pdf  content/auswertung.tex
	content/diskussion.tex content/durchfuehrung.tex content/fehlerrechung.tex
	content/theorie.tex

build/main.pdf: FORCE | build
	  TEXINPUTS="$(call translate,build:)" \
	  BIBINPUTS=build: \
	  max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	main.tex

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
