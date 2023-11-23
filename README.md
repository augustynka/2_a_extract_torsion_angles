# 2_a_extract_torsion_angles
Narzędzie korzysta z bibliotek BioPython Bio.PDB (PDBList, PDBParser, Bio.PDB.vectors - calc_dihedral), pandas, numpy

_Funkcje_ 
**calculate_torsion_angle(atom1,atom2,atom3,atom4)**
funkcja korzysta z współrzędnych wektorowych 4 atomów w celu przeliczenia kątów torsyjnych między płaszczyznami, które te kąty wyznaczają (pierwsza płaszczyzną - atomy 1,2,3- druga płaszczyzną - atomy 2,3,4)
Żeby przeliczyć kąty torsyjne funkcja korzysta z gotowego narzędzia biblioteki (Bio.PDB.vectors.calc_dihedral()), które przyjmuje jako argumenty wartości wektorów wybranych wektorów. Narzędzie to zwraca wartość w radianach, dlatego finalnym krokiem jest przeliczenie wartości na stopnie.
**extract_torsion_angles(structure)**
funkcja tworzy pustą tablicę danych kątów torsyjnych
w serii pętli (w każdym przefiltrowanym łańcuchy w każdym modelu struktury - filtrowanie polega na pominięciu heteroatomów - wybranie reszt, których id zaczyna się od spacji) dla każdej resztynukleotydowej (A,U,G lub C) funkcja tworzy zbiór dla kątów, następnie:
kąt alfa
dla każdej poza pierwszą resztą
pobiera poprzednią resztę, z której potrzebny jest pierwszy atom, następnie według wzoru stosuje funkcję calculate_torsion_angle na atomach: O3(poprzedniej reszty), P(aktualnej), O5 (aktualnej), C5(aktualnej)
kąt beta według atomów P, O5, C5, C4
kąt gamma według atomów O5, C5,C4,C4
kąt delta według atomów C5, C4, C3,O3
kąty dla cukru:
v0 - atomy C4, O4,C1,C2
v1 - atomy O4, C1, C2, C3
v2 - atomy C1, C2, C3, C4
v3 - atomy C2, C3, C4, O4
v4 - atomy C3, C4, O4, C1
kąt chi dla zasady azotowej:
utworzony jest słownik, który mapuje każdą zasady do dwóch atomów, które wchodzą w jej skład
pobiera się atomy z wykorzystaniem słownika i nazwy reszty azotowej, następnie na podstawie atomów O4, C1, [uzyskanego azotu], [uzyskanego węgla] oblicza się kąt torsyjny chi
kąty dodane są do listy torsion_angles_data z wyjątkiem, który sygnalizuje braki w atomach do wyliczania kątów
Całość funkcji zwraca ramkę danych.
**load_structure(pdb_id)**
na podstawie podanego pdb_id pobiera plik .pdb, następnie zwraca strukturę załadowaną do programu za pomocą gotowego parseraPDB

_Działanie programu:_
pdb_id jest ustawione w kodzie [l.87]
wczytywanie struktury
ekstrakcja kątów torsyjnych
zapis ramki danych do pliku .csv w postaci macierzy

