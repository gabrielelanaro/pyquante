AD = ~/Python/asciidoc-7.0.4/asciidoc.py

all: pyquante.html userguide.html

pyquante.html: pyquante.txt
	$(AD) $^

userguide.html: userguide.txt
	$(AD) $^

test:
	open pyquante.html userguide.html

install:
	scp pyquante.html userguide.html rpmuller@pyquante.sf.net:pyquante/htdocs