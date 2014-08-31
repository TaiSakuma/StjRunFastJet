StjRunFastJet
=============

A wrapper library for FastJet used for the STAR experiment

### How to download, compile, and install

Befor installing this librarly, you need to install https://github.com/TaiSakuma/StJetMakerTree and [FastJet 2](http://fastjet.fr/tools-v2.html)

Checkout

    git clone git@github.com:TaiSakuma/StjRunFastJet.git
    cd StjRunFastJet

Run Autotools

    aclocal -I ./m4
    glibtoolize --copy --force
    autoconf
    automake --add-missing --copy --foreign

Configure (replace PREFIX with the path to which you want to install files, e.g., $HOME)

    CPPFLAGS="-I(path to fastjet include files) -I(path to SetJetMakerTree include files" ./configure --prefix=PREFIX

Compile

    make

Install

    make install

this will install files in PREFIX/include and PREFIX/lib
