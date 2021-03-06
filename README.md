# Vigenère Cipher En-/Decryption Program

General reading from `stdin`, capable of:

* Encrypt with a given key
* Decrypt with a given key
* Decrypt by guessing the key via frequency analysis (user `--plot` to generate lots of graphs on the way)


Several tweaking parameters can be given. See `--help` for details.

# Example

Taken from [here](http://www.simonsingh.net/The_Black_Chamber/vigenere_cracking_tool.html):

> RIKVBIYBITHUSEVAZMMLTKASRNHPNPZICSWDSVMBIYFQEZUBZPBRGYNTBURMBECZQKBMBPAWIXSOFNUZECNRAZFPHIYBQEOCTTIOXKUNOHMRGCNDDXZWIRDVDRZYAYYICPUYDHCKXQIECIEWUICJNNACSAZZZGACZHMRGXFTILFNNTSDAFGYWLNICFISEAMRMORPGMJLUSTAAKBFLTIBYXGAVDVXPCTSVVRLJENOWWFINZOWEHOSRMQDGYSDOPVXXGPJNRVILZNAREDUYBTVLIDLMSXKYEYVAKAYBPVTDHMTMGITDZRTIOVWQIECEYBNEDPZWKUNDOZRBAHEGQBXURFGMUECNPAIIYURLRIPTFOYBISEOEDZINAISPBTZMNECRIJUFUCMMUUSANMMVICNRHQJMNHPNCEPUSQDMIVYTSZTRGXSPZUVWNORGQJMYNLILUKCPHDBYLNELPHVKYAYYBYXLERMMPBMHHCQKBMHDKMTDMSSJEVWOPNGCJMYRPYQELCDPOPVPBIEZALKZWTOPRYFARATPBHGLWWMXNHPHXVKBAANAVMNLPHMEMMSZHMTXHTFMQVLILOVVULNIWGVFUCGRZZKAUNADVYXUDDJVKAYUYOWLVBEOZFGTHHSPJNKAYICWITDARZPVU

```
cat above | python vigenere.py --decrypt --plot
```

```Keylength probably 7
Possible passwords:
VIRTUAL
```
![alt text](https://github.com/smoneck/vigenere/blob/master/vigenere_factors.png "Keylength analysis example")

# Pitfalls

1. Just takes CAPITAL LETTERS (since used for [Krypton](http://overthewire.org/wargames/krypton/krypton4.html)). Use e.g. `tr [a-z] [A-Z]` for translation
2. Decryption out of the box just works for *english* material. Adjust monogram frequencies hardcoded in the head of the file for other languages.
