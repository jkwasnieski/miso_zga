"""
Object containing output from Data from miso_vs_miso.miso_bf
JK 180111
"""
import re
import pdb


class MisoEvent(object):
    """
    Object that parses information in miso_vs_miso.miso_bf output file.
    Takes in a single line of file along with file header. 
    Includes methods to return event name and coordinates of each event.
    """

    def __init__(self, line, header):
        if not isinstance(line, str):
            raise ValueError(
                "Must initialize MisoCompare object with a miso_bf file line")
        self._name = None
        self._posterior_mean = None
        self._ci = None
        self._diff = None
        self._bayes = None
        self._isoforms = []
        self._counts = None
        self._assigned_counts = None
        self._chrom = None
        self._strand = None
        self._coords = None
        self.parse_table(line, header)
        self._gene = None

    def __str__(self):
        """ Return Miso Compare object as a string """
        summary = "{event_name}\t{chrom}\t{strand}\t{bayes}\t{diff}\t{isoforms}".format(
            event_name=self._name,
            chrom=self._chrom,
            strand=self._strand,
            bayes=str(self._bayes),
            diff=str(self._diff),
            isoforms="\t".join(self._isoforms))
        return summary + "\t{g}".format(g=self._gene) if self._gene else summary

    def passes_bayes_filter(self, cutoff=100.0):
        """ Returns True if Bayes value is greater than cutoff """
        return True if self._bayes > float(cutoff) else False

    def has_negative_diff(self):
        """ Returns True if diff value is negative """
        return True if self._diff < 0 else False

    def name(self):
        """Returns event name"""
        return self._name

    def gene(self):
        """Returns the trascript's gene name"""
        return self._gene

    def get_event_coords(self):
        """ Returns a list of coordinates that compose the event name """
        return [coord for coord in self.parse_event_name() if coord.isdigit()]

    def get_full_coords(self):
        """ Returns a list of event name coordinates in form chrom:strand:coordinate """
        return [":".join([self._chrom, self._strand, coord]) for coord in self.get_event_coords()]

    def add_gene(self, gene):
        """ Adds gene name attribute """
        self._gene = gene
        return None

    def parse_table(self, line, header):
        """ Create instance of MisoCompare object from miso_vs_miso.miso_bf file """
        # pdb.set_trace()
        values = line.rstrip().split()
        if not len(header) == len(values):
            raise ValueError("{h} number of fields in header but {v} number of fields in "
                             "line".format(h=str(len(header)), v=str(len(values))))
        events = {header[index]: values[index] for index in range(len(header))}
        self._name = events['event_name']
        self._posterior_mean = (events['sample1_posterior_mean'], events[
                                'sample2_posterior_mean'])
        self._ci = ([events['sample1_ci_low'],
                     events['sample1_ci_high']],
                    [events['sample2_ci_low'],
                     events['sample2_ci_high']])
        self._diff = float(events['diff'])
        self._bayes = float(events['bayes_factor'])
        self._isoforms = self.parse_isoforms(events['isoforms'])
        self._counts = (events['sample1_counts'], events['sample2_counts'])
        self._assigned_counts = (events['sample1_assigned_counts'], events[
                                 'sample2_assigned_counts'])
        self._chrom = events['chrom']
        self._strand = events['strand']
        self._coords = (events['mRNA_starts'], events['mRNA_ends'])

    def parse_isoforms(self, isoforms):
        """ Return isoform name parsed into a list of component words """
        return [item for item in isoforms.split("'") if len(item) > 1]

    def parse_event_name(self):
        """ Return event name as a list of chromosome and splice site words
        """
        return re.findall(r"[\w]+", self._name)
