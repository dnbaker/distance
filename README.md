### Distance

This project contains a set of distance metrics for comparing sets, and probably a good starting point for developing new metrics.

#### Sets

Each is templated by float type (to produce the desired floating-point type), and consumes 3 arguments: shared, unique to A, and unique to B.

The classic Jaccard is `a / (a+b+c)`, for instance, while the Dice metric values shared elements over unshared elements.

We support the following from [this paper](http://www.iiisci.org/journal/CV$/sci/pdfs/GS315JG.pdf)\*:

1. Jaccard
2. Dice
3. 3W Jaccard
4. SokalSneath
5. Cosine
6. OchiaiI
7. SorgenFrei
8. MountFord
9. McConnaughey
10. Otsuka
11. KulczynskiII
12. DriverKroeber
13. Johnson
14. Simpson
15. BraunBauquet
16. Hamming
17. Euclid
18. LanceWilliams
19. Hellinger
20. Chord
21. Minkowski

\* Full Citation: `Choi, Seung-Seok, Sung-Hyuk Cha, and Charles C. Tappert. "A survey of binary similarity and distance measures." Journal of Systemics, Cybernetics and Informatics 8.1 (2010): 43-48.

#### Real-valued distances

These aren't yet provided (though blaze-lib, Eigen, libtorch, and many others do) currently, but could be branched into, though we expect to work more with metric learning for real-valued data.
