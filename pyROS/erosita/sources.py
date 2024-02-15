"""
eROSITA source representations.
"""

class eROSITASource:
    """
    Abstract class representation of eROSITA sources from catalogs.
    """
    _standard_keys = []

    def __init__(self,**kwargs):

        # adding provided attributes
        for k,v in kwargs.items():
            setattr(self,k,v)

        for k in self._standard_keys:
            if k not in kwargs:
                setattr(self,k,None)

    @classmethod
    def from_pandas(cls,pandas_df):
        columns = pandas_df.columns

        out = []
        for k in range(len(pandas_df)):
            out.append(cls(**{col:data for col,data in zip(columns,pandas_df.iloc[k,:])}))

        return out

    def __str__(self):
        return f"<eROSITA Source: {self.IAUNAME}>"

    def __repr__(self):
        return self.__str__()





class eRASS1Source(eROSITASource):
    """
    Source corresponding to the eRASS1 data release objects.
    """