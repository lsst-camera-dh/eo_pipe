from collections import defaultdict
import pickle
import pandas as pd


__all__ = ['ETResults']


class ETResults(dict):
    """
    dict subclass to retrieve and provided access to harnessed job
    results from the eT database for a specified run.

    The keys are the schema names and each dict value is a pandas
    dataframe containing the results with a column for each schema
    entry.
    """
    amp_names = 'C10 C11 C12 C13 C14 C15 C16 C17 C07 C06 C05 C04 C03 C02 C01 C00'.split()
    wf_amp_names = 'C10 C11 C12 C13 C14 C15 C16 C17'.split()

    def __init__(self, run, user='ccs', prodServer=True):
        """
        Parameters
        ----------
        run: str
            Run number.  If it ends in 'D', the Dev database will be
            queried.
        user: str ['ccs']
            User id to pass to the etravelerAPI.connection.Connnection
            initializer.
        prodServer: bool [True]
            Flag to use the prod or dev eT server.
        """
        from eTraveler.clientAPI.connection import Connection
        super(ETResults, self).__init__()
        db_name = 'Dev' if run.endswith('D') else 'Prod'
        conn = Connection(user, db_name, prodServer=prodServer)
        self.results = conn.getRunResults(run=run)
        self._extract_schema_values()

    def _extract_schema_values(self):
        steps = self.results['steps']
        for jobname, schema_data in steps.items():
            for schema_name, entries in schema_data.items():
                schema_data = defaultdict(list)
                for entry in entries[1:]:
                    for colname, value in entry.items():
                        schema_data[colname].append(value)
                self[schema_name] = pd.DataFrame(data=schema_data)

    def get_amp_data(self, schema_name, field_name):
        df = self[schema_name]
        amp_data = defaultdict(dict)
        for i in range(len(df)):
            row = df.iloc[i]
            det_name = '_'.join((row.raft, row.slot))
            if 'SW' in det_name:
                amp_data[det_name][self.wf_amp_names[row.amp-1]] \
                    = row[field_name]
            else:
                amp_data[det_name][self.amp_names[row.amp-1]] = row[field_name]
        return amp_data

    def get_amp_gains(self, det_name, schema_name='fe55_BOT_analysis'):
        gain_field = {'fe55_BOT_analysis': 'gain',
                      'ptc_BOT': 'ptc_gain'}
        amp_data = self.get_amp_data(schema_name, gain_field[schema_name])
        return {i+1: amp_data[det_name][amp_name] for i, amp_name
                in enumerate(amp_data[det_name])}

    def to_pickle(self, outfile):
        with open(outfile, 'wb') as fobj:
            pickle.dump(self, fobj)

    @staticmethod
    def read_pickle(pickle_file):
        with open(pickle_file, 'rb') as fobj:
            return pickle.load(fobj)
