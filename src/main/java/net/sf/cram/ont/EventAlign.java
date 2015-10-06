package net.sf.cram.ont;

import java.util.HashMap;
import java.util.Map;

import junit.framework.Assert;

public class EventAlign {
	public enum EA {
		contig, position, reference_kmer, read_index, strand, event_index, event_level_mean, event_length, model_kmer, model_mean, model_stdv, model_name
	}

	public static final byte[] m = new byte[256];
	static {
		m['A'] = 0;
		m['C'] = 1;
		m['G'] = 2;
		m['T'] = 3;
	}
	public static final byte[] b = new byte[4];
	static {
		b[0] = 'A';
		b[1] = 'C';
		b[2] = 'G';
		b[3] = 'T';
	}
	public static final String[] kmers = new String[1024];
	static {
		for (int i = 0; i < 1024; i++) {
			byte[] bases = new byte[5];
			bases[0] = b[(i >> 8) & 3];
			bases[1] = b[(i >> 6) & 3];
			bases[2] = b[(i >> 4) & 3];
			bases[3] = b[(i >> 2) & 3];
			bases[4] = b[i & 3];
			kmers[i] = new String(bases);
		}
	}

	public static final Map<String, Integer> codes = new HashMap<String, Integer>();
	static {
		for (int i = 0; i < 1024; i++) {
			codes.put(kmers[i], i);
		}
	}

	static {
		for (int i = 0; i < 1024; i++) {
			String kmer = kmers[i];
			Assert.assertNotNull(kmer);

			int code = codes.get(kmer);
			Assert.assertTrue(kmer, code >= 0);
			Assert.assertTrue(kmer, code <= 1024);
			Assert.assertEquals(kmer, i, code);
		}
	}

}
