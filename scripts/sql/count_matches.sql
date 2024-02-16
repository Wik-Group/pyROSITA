CREATE TABLE counts_NED
    AS SELECT UID,IAUNAME,RA,DEC,(SELECT COUNT(XREF_NED.UID) FROM XREF_NED WHERE XREF_NED.UID == eROSITA.UID) as N FROM eROSITA GROUP BY UID;
CREATE TABLE counts_SIMBAD
    AS SELECT UID,IAUNAME,RA,DEC,(SELECT COUNT(XREF_SIMBAD.UID) FROM XREF_SIMBAD WHERE XREF_SIMBAD.UID == eROSITA.UID) AS N FROM eROSITA GROUP BY UID;
CREATE TABLE counts
    AS SELECT UID,IAUNAME,RA,DEC, (SELECT N FROM counts_SIMBAD WHERE counts_SIMBAD.UID=UID) + (SELECT N FROM counts_NED WHERE counts_NED.UID=UID) AS N FROM eROSITA GROUP BY UID
