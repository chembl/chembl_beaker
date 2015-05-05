BEAKER_ROOT_URL="https://www.ebi.ac.uk/chembl/api/utils/"

!addHs
curl -X GET ${BEAKER_ROOT_URL}addHs/$(cat addHs.mol | base64 -w 0 | tr "+/" "-_")
curl -X GET ${BEAKER_ROOT_URL}addHs/$(cat addHs.mol | base64 -w 0 | tr "+/" "-_")?addCoords=1
curl -X POST --data-binary @addHs.mol ${BEAKER_ROOT_URL}addHs
curl -X POST -F "file=@addHs.mol" -F "addCoords=1" ${BEAKER_ROOT_URL}addHs

!atomIsInRing
curl -X GET ${BEAKER_ROOT_URL}atomIsInRing/$(cat addHs.mol | base64 -w 0 | tr "+/" "-_")/1/6
curl -X POST -F "file=@addHs.mol" -F "index=1" -F "size=6" ${BEAKER_ROOT_URL}atomIsInRing

!atomRings
curl -X GET ${BEAKER_ROOT_URL}atomRings/$(cat rings.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @rings.mol ${BEAKER_ROOT_URL}atomRings
curl -X POST -F "file=@rings.mol" ${BEAKER_ROOT_URL}atomRings

!bondIsInRing
curl -X GET ${BEAKER_ROOT_URL}bondIsInRing/$(cat rings.mol | base64 -w 0 | tr "+/" "-_")/1/3
curl -X POST -F "file=@rings.mol" -F "index=1" -F "size=3" ${BEAKER_ROOT_URL}bondIsInRing

!bondRings
curl -X GET ${BEAKER_ROOT_URL}bondRings/$(cat rings.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @rings.mol ${BEAKER_ROOT_URL}bondRings
curl -X POST -F "file=@rings.mol" ${BEAKER_ROOT_URL}bondRings

!breakBonds
curl -X GET ${BEAKER_ROOT_URL}breakBonds/$(cat breakBonds.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @breakBonds.mol ${BEAKER_ROOT_URL}breakBonds
curl -X POST -F "file=@breakBonds.mol" ${BEAKER_ROOT_URL}breakBonds

!canonicalizeSmiles
curl -X GET ${BEAKER_ROOT_URL}canonicalizeSmiles/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")
curl -X GET ${BEAKER_ROOT_URL}canonicalizeSmiles/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @aspirin_no_header.smi ${BEAKER_ROOT_URL}canonicalizeSmiles
curl -X POST --data-binary @aspirin_with_header.smi ${BEAKER_ROOT_URL}canonicalizeSmiles
curl -X GET "${BEAKER_ROOT_URL}canonicalizeSmiles/"$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")"?out_delimiter=|&nameHeader=foo"
curl -X POST -F "file=@aspirin_with_header.smi" -F "out_delimiter=|" -F "nameHeader=foo" ${BEAKER_ROOT_URL}canonicalizeSmiles
curl -X POST -F "file=@non_kekule.smi" -F "kekuleSmiles=0" -F "sanitize=0" ${BEAKER_ROOT_URL}canonicalizeSmiles
curl -X GET "${BEAKER_ROOT_URL}canonicalizeSmiles/"$(cat non_kekule.smi | base64 -w 0 | tr "+/" "-_")"?kekuleSmiles=0&sanitize=0"
curl -X POST -F "file=@non_kekule.smi" -F "kekuleSmiles=0" -F "sanitize=1" ${BEAKER_ROOT_URL}canonicalizeSmiles
curl -X GET "${BEAKER_ROOT_URL}canonicalizeSmiles/"$(cat non_kekule.smi | base64 -w 0 | tr "+/" "-_")"?kekuleSmiles=0&sanitize=1"
curl -X POST -F "file=@non_kekule.smi" -F "kekuleSmiles=1" -F "sanitize=1" ${BEAKER_ROOT_URL}canonicalizeSmiles
curl -X GET "${BEAKER_ROOT_URL}canonicalizeSmiles/"$(cat non_kekule.smi | base64 -w 0 | tr "+/" "-_")"?kekuleSmiles=1&sanitize=1"
curl -X POST -F "file=@isomeric.smi" ${BEAKER_ROOT_URL}canonicalizeSmiles
curl -X POST -F "file=@isomeric.smi" -F "isomericSmiles=1" ${BEAKER_ROOT_URL}canonicalizeSmiles
curl -X GET "${BEAKER_ROOT_URL}canonicalizeSmiles/"$(cat isomeric.smi | base64 -w 0 | tr "+/" "-_")"?isomericSmiles=1"

!ctab23D
curl -X GET ${BEAKER_ROOT_URL}ctab23D/$(cat no_coords.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @no_coords.mol ${BEAKER_ROOT_URL}ctab23D
curl -X POST -F "file=@no_coords.mol" ${BEAKER_ROOT_URL}ctab23D

!ctab2image
curl -X GET ${BEAKER_ROOT_URL}ctab2image/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_") > aspirin.png
curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}ctab2image > aspirin.png
curl -X GET ${BEAKER_ROOT_URL}ctab2image/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?computeCoords=0 > aspirin.png
curl -X POST -F "file=@aspirin.mol" -F "computeCoords=0" ${BEAKER_ROOT_URL}ctab2image > aspirin.png
curl -X GET ${BEAKER_ROOT_URL}ctab2image/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?atomMapNumber=1 > aspirin.png
curl -X POST -F "file=@aspirin.mol" -F "atomMapNumber=1" ${BEAKER_ROOT_URL}ctab2image > aspirin.png
curl -X GET ${BEAKER_ROOT_URL}ctab2image/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?legend=aspirin > aspirin.png
curl -X POST -F "file=@aspirin.mol" -F "legend=aspirin" ${BEAKER_ROOT_URL}ctab2image > aspirin.png
curl -X POST -F "file=@mcs.sdf" -F "legend=foo|bar|bla" ${BEAKER_ROOT_URL}ctab2image > out.png
curl -X GET "${BEAKER_ROOT_URL}ctab2image/"$(cat mcs.sdf | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla" > out.png
curl -X POST -F "file=@mcs.sdf" -F "legend=foo|bar|bla" -F "computeCoords=0" ${BEAKER_ROOT_URL}ctab2image > out.png
curl -X GET "${BEAKER_ROOT_URL}ctab2image/"$(cat mcs.sdf | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla&computeCoords=0" > out.png
curl -X POST -F "file=@mcs_no_coords.sdf" -F "legend=foo|bar|bla" ${BEAKER_ROOT_URL}ctab2image > out.png
curl -X GET "${BEAKER_ROOT_URL}ctab2image/"$(cat mcs.sdf | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla" > out.png
curl -X GET ${BEAKER_ROOT_URL}ctab2image/$(cat mcs_no_coords.sdf | base64 -w 0 | tr "+/" "-_")?legend=foo > out.png
curl -X POST -F "file=@mcs.sdf" -F "legend=foo" ${BEAKER_ROOT_URL}ctab2image > out.png
curl -X POST -F "file=@mcs.sdf" -F "legend=foo|bar|bla" -F "size=400" ${BEAKER_ROOT_URL}ctab2image > out.png
curl -X GET ${BEAKER_ROOT_URL}ctab2image/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?size=400 > aspirin.png
curl -X POST -F "file=@aspirin.mol" -F "size=400" ${BEAKER_ROOT_URL}ctab2image > aspirin.png

!ctab2inchi
curl -X GET ${BEAKER_ROOT_URL}ctab2inchi/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}ctab2inchi
curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}ctab2inchi

!ctab2inchiKey
curl -X GET ${BEAKER_ROOT_URL}ctab2inchiKey/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}ctab2inchiKey
curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}ctab2inchiKey

!ctab2json
!

!ctab2smiles
curl -X GET ${BEAKER_ROOT_URL}ctab2smiles/$(cat isomeric.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST -F "file=@isomeric.mol" ${BEAKER_ROOT_URL}ctab2smiles
curl -X GET ${BEAKER_ROOT_URL}ctab2smiles/$(cat isomeric.mol | base64 -w 0 | tr "+/" "-_")?isomericSmiles=1
curl -X POST -F "file=@isomeric.mol" -F "isomericSmiles=1" ${BEAKER_ROOT_URL}ctab2smiles
curl -X GET "${BEAKER_ROOT_URL}ctab2smiles/"$(cat non_kekule.mol | base64 -w 0 | tr "+/" "-_")"?kekuleSmiles=0&sanitize=1"
curl -X POST -F "file=@non_kekule.mol" -F "kekuleSmiles=0" -F "sanitize=1" ${BEAKER_ROOT_URL}ctab2smiles
curl -X GET "${BEAKER_ROOT_URL}ctab2smiles/"$(cat non_kekule.mol | base64 -w 0 | tr "+/" "-_")"?kekuleSmiles=0&sanitize=0"
curl -X POST -F "file=@non_kekule.mol" -F "kekuleSmiles=0" -F "sanitize=0" ${BEAKER_ROOT_URL}ctab2smiles
curl -X GET "${BEAKER_ROOT_URL}ctab2smiles/"$(cat non_kekule.mol | base64 -w 0 | tr "+/" "-_")"?kekuleSmiles=1&sanitize=1"
curl -X POST -F "file=@non_kekule.mol" -F "kekuleSmiles=1" -F "sanitize=1" ${BEAKER_ROOT_URL}ctab2smiles
curl -X GET "${BEAKER_ROOT_URL}ctab2smiles/"$(cat explicitHs.mol | base64 -w 0 | tr "+/" "-_")"?removeHs=0"
curl -X POST -F "file=@explicitHs.mol" -F "removeHs=0" ${BEAKER_ROOT_URL}ctab2smiles

!ctab2svg
curl -X GET ${BEAKER_ROOT_URL}ctab2svg/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_") > aspirin.svg
curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}ctab2svg > aspirin.svg
curl -X GET ${BEAKER_ROOT_URL}ctab2svg/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?computeCoords=0 > aspirin.svg
curl -X POST -F "file=@aspirin.mol" -F "computeCoords=0" ${BEAKER_ROOT_URL}ctab2svg > aspirin.svg
curl -X GET ${BEAKER_ROOT_URL}ctab2svg/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?atomMapNumber=1 > aspirin.svg
curl -X POST -F "file=@aspirin.mol" -F "atomMapNumber=1" ${BEAKER_ROOT_URL}ctab2svg > aspirin.svg
curl -X GET ${BEAKER_ROOT_URL}ctab2svg/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?legend=aspirin > aspirin.svg
curl -X POST -F "file=@aspirin.mol" -F "legend=aspirin" ${BEAKER_ROOT_URL}ctab2svg > aspirin.svg
curl -X GET ${BEAKER_ROOT_URL}ctab2svg/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")?size=400 > aspirin.svg
curl -X POST -F "file=@aspirin.mol" -F "size=400" ${BEAKER_ROOT_URL}ctab2svg > aspirin.svg

!descriptors
curl -X GET ${BEAKER_ROOT_URL}descriptors/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}descriptors
curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}descriptors

!getNumAtoms
curl -X GET ${BEAKER_ROOT_URL}getNumAtoms/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}getNumAtoms
curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}getNumAtoms

!getNumBonds
curl -X GET ${BEAKER_ROOT_URL}getNumBonds/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}getNumBonds
curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}getNumBonds

!hydrogenize
curl -X POST -H "Content-Type: application/json; charset=UTF-8" -d @metane.json ${BEAKER_ROOT_URL}hydrogenize

!image2ctab
curl -X POST --data-binary @mol.jpg ${BEAKER_ROOT_URL}image2ctab
curl -X POST -F "file=@mol.jpg" ${BEAKER_ROOT_URL}image2ctab
curl -X POST --data-binary @mol.png ${BEAKER_ROOT_URL}image2ctab
curl -X POST -F "file=@mol.png" ${BEAKER_ROOT_URL}image2ctab
curl -X GET ${BEAKER_ROOT_URL}image2ctab/$(cat mol.jpg | base64 -w 0 | tr "+/" "-_")
curl -X GET ${BEAKER_ROOT_URL}image2ctab/$(cat mol.png | base64 -w 0 | tr "+/" "-_")

!image2smiles
curl -X POST --data-binary @mol.jpg ${BEAKER_ROOT_URL}image2smiles
curl -X POST -F "file=@mol.jpg" ${BEAKER_ROOT_URL}image2smiles
curl -X POST --data-binary @mol.png ${BEAKER_ROOT_URL}image2smiles
curl -X POST -F "file=@mol.png" ${BEAKER_ROOT_URL}image2smiles
curl -X GET ${BEAKER_ROOT_URL}image2smiles/$(cat mol.jpg | base64 -w 0 | tr "+/" "-_")
curl -X GET ${BEAKER_ROOT_URL}image2smiles/$(cat mol.png | base64 -w 0 | tr "+/" "-_")

!inchi2ctab
curl -X GET ${BEAKER_ROOT_URL}inchi2ctab/$(cat aspirin.inchi | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @aspirin.inchi ${BEAKER_ROOT_URL}inchi2ctab
curl -X POST -F "file=@aspirin.inchi" ${BEAKER_ROOT_URL}inchi2ctab

!inchi2inchiKey
curl -X GET ${BEAKER_ROOT_URL}inchi2inchiKey/$(cat aspirin.inchi | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @aspirin.inchi ${BEAKER_ROOT_URL}inchi2inchiKey
curl -X POST -F "file=@aspirin.inchi" ${BEAKER_ROOT_URL}inchi2inchiKey

!inchi2svg
curl -X GET ${BEAKER_ROOT_URL}inchi2svg/$(cat aspirin.inchi | base64 -w 0 | tr "+/" "-_") > aspirin.svg
curl -X POST --data-binary @aspirin.inchi ${BEAKER_ROOT_URL}inchi2svg > aspirin.svg
curl -X GET ${BEAKER_ROOT_URL}inchi2svg/$(cat aspirin.inchi | base64 -w 0 | tr "+/" "-_")?legend=aspirin > aspirin.svg
curl -X POST -F "file=@aspirin.inchi" -F "legend=aspirin" ${BEAKER_ROOT_URL}inchi2svg > aspirin.svg
curl -X GET ${BEAKER_ROOT_URL}inchi2svg/$(cat aspirin.inchi | base64 -w 0 | tr "+/" "-_")?size=400 > aspirin.svg
curl -X POST -F "file=@aspirin.inchi" -F "size=400" ${BEAKER_ROOT_URL}inchi2svg > aspirin.svg

!kekulize
curl -X GET ${BEAKER_ROOT_URL}kekulize/$(cat aromatic.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @aromatic.mol ${BEAKER_ROOT_URL}kekulize
curl -X POST -F "file=@aromatic.mol" ${BEAKER_ROOT_URL}kekulize

!logP
curl -X GET ${BEAKER_ROOT_URL}logP/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}logP
curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}logP

!mcs
curl -X GET ${BEAKER_ROOT_URL}mcs/$(cat mcs.sdf | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @mcs.sdf ${BEAKER_ROOT_URL}mcs
curl -X POST -F "file=@mcs.sdf" ${BEAKER_ROOT_URL}mcs

!molExport
curl -X POST -H "Content-Type: application/json; charset=UTF-8" -d @molExport.json ${BEAKER_ROOT_URL}molExport

!molWt
curl -X GET ${BEAKER_ROOT_URL}molWt/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}molWt
curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}molWt

!neutralise
curl -X GET ${BEAKER_ROOT_URL}neutralise/$(cat neutralise.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @neutralise.mol ${BEAKER_ROOT_URL}neutralise
curl -X POST -F "file=@neutralise.mol" ${BEAKER_ROOT_URL}neutralise

!numAtomRings
curl -X GET ${BEAKER_ROOT_URL}numAtomRings/$(cat rings.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @rings.mol ${BEAKER_ROOT_URL}numAtomRings
curl -X POST -F "file=@rings.mol" ${BEAKER_ROOT_URL}numAtomRings

!numBondRings
curl -X GET ${BEAKER_ROOT_URL}numBondRings/$(cat rings.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @rings.mol ${BEAKER_ROOT_URL}numBondRings
curl -X POST -F "file=@rings.mol" ${BEAKER_ROOT_URL}numBondRings

!numRings
curl -X GET ${BEAKER_ROOT_URL}numRings/$(cat rings.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @rings.mol ${BEAKER_ROOT_URL}numRings
curl -X POST -F "file=@rings.mol" ${BEAKER_ROOT_URL}numRings

!removeHs - przydalby sie przyklad z implicitOnly, w ktorym jednak cos zostaje usuniete...
curl -X GET ${BEAKER_ROOT_URL}removeHs/$(cat removeHs.mol | base64 -w 0 | tr "+/" "-_")
curl -X GET ${BEAKER_ROOT_URL}removeHs/$(cat removeHs.mol | base64 -w 0 | tr "+/" "-_")?implicitOnly=1
curl -X POST --data-binary @removeHs.mol ${BEAKER_ROOT_URL}removeHs
curl -X POST -F "file=@removeHs.mol" -F "implicitOnly=1" ${BEAKER_ROOT_URL}removeHs

!rules
curl -X GET ${BEAKER_ROOT_URL}rules/$(cat rules.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @rules.mol ${BEAKER_ROOT_URL}rules
curl -X POST -F "file=@rules.mol" ${BEAKER_ROOT_URL}rules

!sanitize
curl -X GET ${BEAKER_ROOT_URL}sanitize/$(cat aromatic.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @aromatic.mol ${BEAKER_ROOT_URL}sanitize
curl -X POST -F "file=@aromatic.mol" ${BEAKER_ROOT_URL}sanitize

!sdf2SimilarityMap
curl -X GET ${BEAKER_ROOT_URL}sdf2SimilarityMap/$(cat sim.sdf | base64 -w 0 | tr "+/" "-_") > sim.png
curl -X GET "${BEAKER_ROOT_URL}sdf2SimilarityMap/"$(cat sim.sdf | base64 -w 0 | tr "+/" "-_")"?width=500&height=500" > sim.png
curl -X GET "${BEAKER_ROOT_URL}sdf2SimilarityMap/"$(cat sim.sdf | base64 -w 0 | tr "+/" "-_")"?width=500&height=500&fingerprint=tt" > sim.png
curl -X GET "${BEAKER_ROOT_URL}sdf2SimilarityMap/"$(cat sim.sdf | base64 -w 0 | tr "+/" "-_")"?width=500&height=500&fingerprint=ap" > sim.png
curl -X POST --data-binary @sim.sdf ${BEAKER_ROOT_URL}sdf2SimilarityMap > sim.png
curl -X POST -F "file=@sim.sdf" -F "width=500" -F "height=500" -F "fingerprint=tt" ${BEAKER_ROOT_URL}sdf2SimilarityMap > sim.png
curl -X POST -F "file=@sim.sdf" -F "width=500" -F "height=500" -F "fingerprint=ap" ${BEAKER_ROOT_URL}sdf2SimilarityMap > sim.png

!smiles2SimilarityMap
curl -X GET ${BEAKER_ROOT_URL}smiles2SimilarityMap/$(cat sim.smi | base64 -w 0 | tr "+/" "-_") > sim.png
curl -X GET "${BEAKER_ROOT_URL}smiles2SimilarityMap/"$(cat sim.smi | base64 -w 0 | tr "+/" "-_")"?width=500&height=500" > sim.png
curl -X GET "${BEAKER_ROOT_URL}smiles2SimilarityMap/"$(cat sim.smi | base64 -w 0 | tr "+/" "-_")"?width=500&height=500&fingerprint=tt" > sim.png
curl -X GET "${BEAKER_ROOT_URL}smiles2SimilarityMap/"$(cat sim.smi | base64 -w 0 | tr "+/" "-_")"?width=500&height=500&fingerprint=ap" > sim.png
curl -X POST --data-binary @sim.smi ${BEAKER_ROOT_URL}smiles2SimilarityMap > sim.png
curl -X POST -F "file=@sim.smi" -F "width=500" -F "height=500" -F "fingerprint=tt" ${BEAKER_ROOT_URL}smiles2SimilarityMap > sim.png
curl -X POST -F "file=@sim.smi" -F "width=500" -F "height=500" -F "fingerprint=ap" ${BEAKER_ROOT_URL}smiles2SimilarityMap > sim.png

!sdf2fps
curl -X GET ${BEAKER_ROOT_URL}sdf2fps/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
curl -X GET "${BEAKER_ROOT_URL}sdf2fps/"$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")"?n_bits=1024&radius=5"
curl -X GET "${BEAKER_ROOT_URL}sdf2fps/"$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")"?type=pair"
curl -X GET "${BEAKER_ROOT_URL}sdf2fps/"$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")"?type=maccs"
curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}sdf2fps
curl -X POST -F "file=@aspirin.mol" -F "n_bits=1024" -F "radius=5" ${BEAKER_ROOT_URL}sdf2fps
curl -X POST -F "file=@aspirin.mol" -F "type=maccs" ${BEAKER_ROOT_URL}sdf2fps

!smiles23D
curl -X GET ${BEAKER_ROOT_URL}smiles23D/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @aspirin_with_header.smi ${BEAKER_ROOT_URL}smiles23D
curl -X POST -F "file=@aspirin_with_header.smi" ${BEAKER_ROOT_URL}smiles23D
curl -X GET ${BEAKER_ROOT_URL}smiles23D/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @aspirin_no_header.smi ${BEAKER_ROOT_URL}smiles23D
curl -X POST -F "file=@aspirin_no_header.smi" ${BEAKER_ROOT_URL}smiles23D

!smiles2ctab
curl -X GET ${BEAKER_ROOT_URL}smiles2ctab/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")
curl -X POST -F "file=@aspirin_with_header.smi" ${BEAKER_ROOT_URL}smiles2ctab
curl -X GET ${BEAKER_ROOT_URL}smiles2ctab/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")
curl -X POST -F "file=@aspirin_no_header.smi" ${BEAKER_ROOT_URL}smiles2ctab
curl -X GET "${BEAKER_ROOT_URL}smiles2ctab/"$(cat rules.smi | base64 -w 0 | tr "+/" "-_")"?computeCoords=0"
curl -X POST -F "file=@rules.smi" -F "computeCoords=0"  ${BEAKER_ROOT_URL}smiles2ctab
curl -X GET ${BEAKER_ROOT_URL}smiles2ctab/$(cat mcs.smi | base64 -w 0 | tr "+/" "-_")
curl -X POST -F "file=@mcs.smi" ${BEAKER_ROOT_URL}smiles2ctab
curl -X GET ${BEAKER_ROOT_URL}smiles2ctab/$(cat mcs_no_header.smi | base64 -w 0 | tr "+/" "-_")
curl -X POST -F "file=@mcs_no_header.smi" ${BEAKER_ROOT_URL}smiles2ctab

!smiles2image
curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_") > aspirin.png
curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_") > aspirin.png
curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?atomMapNumber=1 > aspirin.png
curl -X POST -F "file=@aspirin_no_header.smi" -F "atomMapNumber=1" ${BEAKER_ROOT_URL}smiles2image > aspirin.png
curl -X POST --data-binary @aspirin_with_header.smi ${BEAKER_ROOT_URL}smiles2image > aspirin.png
curl -X POST --data-binary @aspirin_no_header.smi ${BEAKER_ROOT_URL}smiles2image > aspirin.png
curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")?legend=aspirin > aspirin.png
curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?legend=aspirin > aspirin.png
curl -X POST -F "file=@aspirin_with_header.smi" -F "legend=aspirin" ${BEAKER_ROOT_URL}smiles2image > aspirin.png
curl -X POST -F "file=@aspirin_no_header.smi" -F "legend=aspirin" ${BEAKER_ROOT_URL}smiles2image > aspirin.png
curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")?size=400 > aspirin.png
curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?size=400 > aspirin.png
curl -X POST -F "file=@aspirin_with_header.smi" -F "size=400" ${BEAKER_ROOT_URL}smiles2image > aspirin.png
curl -X POST -F "file=@aspirin_no_header.smi" -F "size=400" ${BEAKER_ROOT_URL}smiles2image > aspirin.png
curl -X POST -F "file=@mcs.smi" -F "legend=foo|bar|bla" ${BEAKER_ROOT_URL}smiles2image > out.png
curl -X POST -F "file=@mcs_no_header.smi" -F "legend=foo|bar|bla" ${BEAKER_ROOT_URL}smiles2image > out.png
curl -X GET "${BEAKER_ROOT_URL}smiles2image/"$(cat mcs.smi | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla" > out.png
curl -X GET "${BEAKER_ROOT_URL}smiles2image/"$(cat mcs_no_header.smi | base64 -w 0 | tr "+/" "-_")"?legend=foo|bar|bla" > out.png
curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat mcs.smi | base64 -w 0 | tr "+/" "-_")?legend=foo > out.png
curl -X GET ${BEAKER_ROOT_URL}smiles2image/$(cat mcs_no_header.smi | base64 -w 0 | tr "+/" "-_")?legend=foo > out.png
curl -X POST -F "file=@mcs.smi" -F "legend=foo" ${BEAKER_ROOT_URL}smiles2image > out.png
curl -X POST -F "file=@mcs_no_header.smi" -F "legend=foo" ${BEAKER_ROOT_URL}smiles2image > out.png
curl -X POST -F "file=@mcs.smi" -F "legend=foo|bar|bla" -F "size=400" ${BEAKER_ROOT_URL}smiles2image > out.png
curl -X POST -F "file=@mcs_no_header.smi" -F "legend=foo|bar|bla" -F "size=400" ${BEAKER_ROOT_URL}smiles2image > out.png

!smiles2inchi
curl -X GET ${BEAKER_ROOT_URL}smiles2inchi/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")
curl -X POST -F "file=@aspirin_with_header.smi" ${BEAKER_ROOT_URL}smiles2inchi
curl -X GET ${BEAKER_ROOT_URL}smiles2inchi/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")
curl -X POST -F "file=@aspirin_no_header.smi" ${BEAKER_ROOT_URL}smiles2inchi
curl -X GET ${BEAKER_ROOT_URL}smiles2inchi/$(cat mcs.smi | base64 -w 0 | tr "+/" "-_")
curl -X POST -F "file=@mcs.smi" ${BEAKER_ROOT_URL}smiles2inchi
curl -X GET ${BEAKER_ROOT_URL}smiles2inchi/$(cat mcs_no_header.smi | base64 -w 0 | tr "+/" "-_")
curl -X POST -F "file=@mcs_no_header.smi" ${BEAKER_ROOT_URL}smiles2inchi

!smiles2inchiKey
curl -X GET ${BEAKER_ROOT_URL}smiles2inchiKey/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")
curl -X POST -F "file=@aspirin_with_header.smi" ${BEAKER_ROOT_URL}smiles2inchiKey
curl -X GET ${BEAKER_ROOT_URL}smiles2inchiKey/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")
curl -X POST -F "file=@aspirin_no_header.smi" ${BEAKER_ROOT_URL}smiles2inchi
curl -X GET ${BEAKER_ROOT_URL}smiles2inchiKey/$(cat mcs.smi | base64 -w 0 | tr "+/" "-_")
curl -X POST -F "file=@mcs.smi" ${BEAKER_ROOT_URL}smiles2inchiKey
curl -X GET ${BEAKER_ROOT_URL}smiles2inchiKey/$(cat mcs_no_header.smi | base64 -w 0 | tr "+/" "-_")
curl -X POST -F "file=@mcs_no_header.smi" ${BEAKER_ROOT_URL}smiles2inchiKey

!smiles2json
!

!smiles2svg
curl -X GET ${BEAKER_ROOT_URL}smiles2svg/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_") > aspirin.svg
curl -X POST --data-binary @aspirin_no_header.smi ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
curl -X GET ${BEAKER_ROOT_URL}smiles2svg/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?atomMapNumber=1 > aspirin.svg
curl -X POST -F "file=@aspirin_no_header.smi" -F "atomMapNumber=1" ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
curl -X GET ${BEAKER_ROOT_URL}smiles2svg/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?legend=aspirin > aspirin.svg
curl -X POST -F "file=@aspirin_no_header.smi" -F "legend=aspirin" ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
curl -X GET ${BEAKER_ROOT_URL}smiles2svg/$(cat aspirin_no_header.smi | base64 -w 0 | tr "+/" "-_")?size=400 > aspirin.svg
curl -X POST -F "file=@aspirin_no_header.smi" -F "size=400" ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
curl -X GET ${BEAKER_ROOT_URL}smiles2svg/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_") > aspirin.svg
curl -X POST --data-binary @aspirin_with_header.smi ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
curl -X GET ${BEAKER_ROOT_URL}smiles2svg/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")?legend=aspirin > aspirin.svg
curl -X POST -F "file=@aspirin_with_header.smi" -F "legend=aspirin" ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg
curl -X GET ${BEAKER_ROOT_URL}smiles2svg/$(cat aspirin_with_header.smi | base64 -w 0 | tr "+/" "-_")?size=400 > aspirin.svg
curl -X POST -F "file=@aspirin_with_header.smi" -F "size=400" ${BEAKER_ROOT_URL}smiles2svg > aspirin.svg

!sssr
curl -X GET ${BEAKER_ROOT_URL}sssr/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}sssr
curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}sssr

!standardise
curl -X GET ${BEAKER_ROOT_URL}standardise/$(cat standardise.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @standardise.mol ${BEAKER_ROOT_URL}standardise
curl -X POST -F "file=@standardise.mol" ${BEAKER_ROOT_URL}standardise

!symmSSSR
curl -X GET ${BEAKER_ROOT_URL}symmSSSR/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}symmSSSR
curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}symmSSSR

!tpsa
curl -X GET ${BEAKER_ROOT_URL}tpsa/$(cat aspirin.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @aspirin.mol ${BEAKER_ROOT_URL}tpsa
curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}tpsa

!unsalt
curl -X GET ${BEAKER_ROOT_URL}unsalt/$(cat unsalt.mol | base64 -w 0 | tr "+/" "-_")
curl -X POST --data-binary @unsalt.mol ${BEAKER_ROOT_URL}unsalt
curl -X POST -F "file=@aspirin.mol" ${BEAKER_ROOT_URL}unsalt

