# Analysis of human ENCODE4 long read data

* [ENCODE cart of human cell line data](https://www.encodeproject.org/carts/723c6f14-e68b-4480-8a61-704a15ac5c7a/)
* [ENCODE cart of human tissue data](https://www.encodeproject.org/carts/26bd2879-329d-4168-98b9-6d132a1aad0f/)
* [ENCODE cart of all human data](https://www.encodeproject.org/carts/829d339c-913c-4773-8001-80130796a367/)

## Download processed files (post-TranscriptClean)
```bash
xargs -L 1 curl -O -J -L < files.txt
```
