"""
test_phylovector
================

Tests for the `phylovector` package.
"""

import pytest
import ete3

# Import the library being tested
import phylovector

# Random trees generated with `ngesh`


@pytest.mark.parametrize(
    "newick",
    [
        "((Kotu:0.0871121,(Ao:0.589289,(Vigep:0.502903,Zupkopap:0.352368):0.0863852):0.0603562):1.24338,(((Diggi:0.52503,Kota:0.264633):0.384149,(Sabafe:0.971275,Tuba:0.971275):0.120429):0.381786,((Kitedo:0.0192851,Vuvu:0.0192851):1.20977,(Marugi:0.703817,Teguo:0.310383):0.525241):0.244433):0.419535);",
        "(Torduo:1.93864,(Veita:0.345091,((Kuodek:0.0170481,Ousu:0.0170481):0.316698,(Rohe:0.217531,Wigisi:0.217531):0.116215):0.504362):1.10053);",
        "((Metazo:0.60289,(Aga:0.261902,(Haholef:0.0836789,Nosula:0.0836789):0.178223):0.340987):1.24527,((Kino:0.245433,Ule:0.323594):0.249263,(Kozo:1.22807,(Dobomto:0.844976,Mu:0.844976):0.383097):0.232809):0.387277);",
        "(Adutam:0.0425664,((Ufi:1.47105,Wizof:1.47105):0.34336,(Ueda:0.0790154,(Tibo:0.532147,((Bavlaz:0.393377,Kumou:0.136403):0.00517799,(Lemkoar:0.39576,Tiki:0.39576):0.00279468):0.133592):0.0348807):1.24738):0.143984);",
        "((Kunih:0.202793,(Lupikzi:0.396125,Nozi:0.0138682):0.13849):0.691951,((Sezaho:0.452701,Zigfei:0.452701):1.18883,(Sudo:0.0704571,(Lono:1.00289,Nerez:0.027156):0.0702991):0.568335):0.21791);",
        "((Huela:0.114706,Nire:0.114706):1.59165,(Sepeo:1.50903,(Dia:0.230105,(Oi:0.0523563,Pugu:0.0491651):0.177748):1.27892):0.197329);",
        "(((Aua:0.365422,Vipe:1.06091):0.441997,((Gonbidri:0.15318,Sifi:0.15318):0.41219,(Rubrai:0.247082,Tedeki:0.247082):0.318288):0.93754):0.393577,(Modo:0.149479,((Ofu:0.131608,Wio:0.369087):0.824902,((Helezum:0.183952,Moze:0.547124):0.200082,(Ubavum:0.564263,Veteltah:0.564263):0.182943):0.446782):0.503512):0.198986);",
        "(Pilo:1.59205,((Uhi:0.216126,(Molkuo:0.201413,(Erie:0.0828869,Oai:0.0828869):0.118527):0.0147131):0.883773,(Zupru:1.00155,(Guna:0.343031,(Nudu:0.836214,(Rozo:0.219106,(Hemif:0.405284,Kegbimlon:0.405284):0.205327):0.225603):0.0435484):0.121785):0.0983528):0.832601);",
        "(Vogab:1.04313,((Ue:0.122688,Vinedfap:0.276122):0.541514,((Bia:0.485216,Maue:0.485216):0.342772,(Pevu:0.0733257,Tutano:0.0733257):0.754663):0.607921):0.346808);",
        "(Lolol:0.900328,(Zua:1.71472,((Keme:1.3712,(Bema:0.080973,Tugugum:0.080973):1.39579):0.071799,(Eki:1.23722,((Aki:0.217514,Desu:0.217514):0.961745,(Zovle:0.819957,(Ezape:0.657111,(Kupge:0.0879704,Tazvi:0.0879704):0.56914):0.162846):0.359303):0.057965):0.311337):0.166161):0.224726);",
    ],
)
def test_from_newick(newick):
    #    print("...", [newick])
    tree_a = ete3.Tree(newick)
    leaves = sorted([leaf.name for leaf in tree_a.iter_leaves()])
    vector = phylovector.tree2vector(tree_a)
    tree_b = phylovector.vector2tree(vector, leaves)
    test_newick = phylovector.sorted_newick(tree_b.write(format=1))

    assert test_newick == newick
