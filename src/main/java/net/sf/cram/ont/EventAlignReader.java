package net.sf.cram.ont;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;
import java.util.zip.GZIPInputStream;

import net.sf.cram.ont.EventAlign.EA;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.SequenceUtil;

public class EventAlignReader {
	private Scanner scanner;
	private Event firstEvent;

	public EventAlignReader(File file) throws IOException {
		super();
		if (file.getName().endsWith(".gz"))
			scanner = new Scanner(new GZIPInputStream(new FileInputStream(file)));
		else
			scanner = new Scanner(file);
		// skip header line:
		if (scanner.hasNextLine()) {
			scanner.nextLine();
			if (scanner.hasNextLine())
				firstEvent = Event.fromString(scanner.nextLine());
		}
	}

	public List<Event> nextRead() {
		if (firstEvent == null)
			return null;

		List<Event> events = new ArrayList<EventAlignReader.Event>();
		events.add(firstEvent);

		firstEvent = null;
		while (scanner.hasNextLine()) {
			String line = scanner.nextLine();
			firstEvent = Event.fromString(line);
			if (!EventAlign.codes.containsKey(firstEvent.model_kmer)) {
				continue;
			}
			// if (firstEvent.strand
			// && !firstEvent.reference_kmer.equals(firstEvent.model_kmer))
			// System.err.println(line);
			Event lastEvent = events.get(events.size() - 1);
			if (firstEvent.read_index != lastEvent.read_index || firstEvent.strand != lastEvent.strand)
				break;
			events.add(firstEvent);
			firstEvent = null;
		}
		if (events.isEmpty())
			return null;

		return events;
	}

	public EASeries transpose(List<Event> events) {
		EASeries series = new EASeries(events.size());
		Event first = events.get(0);
		series.event_index_start = first.event_index;
		series.read_index = first.read_index;
		series.template = first.strand;
		for (int i = 0; i < events.size(); i++) {
			Event event = events.get(i);
			series.event_length[i] = event.event_length;
			series.event_level_mean[i] = event.event_level_mean;
			series.position[i] = event.position;

			String model_kmer = event.strand ? event.model_kmer : SequenceUtil.reverseComplement(event.model_kmer);
			series.model_kmer[i] = EventAlign.codes.get(model_kmer);

			series.scores[i] = score(event.event_level_mean, event.model_mean, event.model_stdv);
		}
		return series;
	}

	public byte[] calculateBases(List<Event> events) {
		byte[] bases = new byte[events.size()];
		int prev = 0;
		List<String> mis = new ArrayList<String>();
		// bases.add(kmers[s.model_kmer[0]].getBytes()[0]);
		// bases.add(kmers[s.model_kmer[0]].getBytes()[1]);
		// bases.add(kmers[s.model_kmer[0]].getBytes()[2]);
		// bases.add(kmers[s.model_kmer[0]].getBytes()[3]);

		Event event = events.get(0);
		int prevPos = event.position;
		String kmer = event.strand ? event.model_kmer : SequenceUtil.reverseComplement(event.model_kmer);
		bases[0] = kmer.getBytes()[4];
		for (int i = 1; i < bases.length; i++) {
			event = events.get(i);
			kmer = event.strand ? event.model_kmer : SequenceUtil.reverseComplement(event.model_kmer);
			byte[] kmerBases = kmer.getBytes();
			int delta = event.position - prevPos;
			switch (delta) {
			case 0:
				bases[i] = kmerBases[4];
				break;
			case 1:
				bases[i] = kmerBases[4];
				break;

			default:
				bases[i] = kmerBases[4];
				break;
			}

			prevPos = event.position;
		}

		return bases;
	}

	public byte[] calculateScores(List<Event> events) {
		byte[] scores = new byte[events.size()];

		int i = 0;
		for (Event event : events)
			scores[i++] = score(event.event_level_mean, event.model_mean, event.model_stdv);

		return scores;
	}

	public String calculateED(List<Event> events) {
		List<String> mis = new ArrayList<String>();
		Event event = events.get(0);
		int prevPos = event.position;
		String kmer = event.strand ? event.model_kmer : SequenceUtil.reverseComplement(event.model_kmer);
		int match = 1;
		for (int i = 1; i < events.size(); i++) {
			event = events.get(i);
			kmer = event.strand ? event.model_kmer : SequenceUtil.reverseComplement(event.model_kmer);

			int delta = event.position - prevPos;
			switch (delta) {
			case 0:
				break;
			case 1:
				match++;
				break;

			default:
				// this captures only differing bases into mis:
				// byte[] em = new byte[delta - 1];
				// for (int j = 0; j < em.length; j++) {
				// em[j] = kmerBases[4 - delta + j + 1];
				// }
				// if (match > 0) {
				// mis.add(match + "M");
				// match = 0;
				// }
				// mis.add("^" + new String(em));

				// the following passage captures the kmer if it is different
				// into mis:
				if (match > 0) {
					mis.add(match + "");
					match = 0;
				}
				mis.add("^" + kmer);
				break;
			}

			prevPos = event.position;
			if (!event.reference_kmer.equals(kmer)) {
				if (match > 0) {
					mis.add(match + "");
					match = 0;
				}
				mis.add(kmer);
			}
		}
		if (match > 0) {
			mis.add(match + "");
			match = 0;
		}

		StringBuffer sb = new StringBuffer();
		for (String m : mis)
			sb.append(m);
		return sb.toString();
	}

	public String calculateMD(List<Event> events) {
		List<String> mis = new ArrayList<String>();
		Event event = events.get(0);
		int prevPos = event.position;
		String kmer = event.strand ? event.model_kmer : SequenceUtil.reverseComplement(event.model_kmer);
		int match = 1;
		for (int i = 1; i < events.size(); i++) {
			event = events.get(i);
			kmer = event.strand ? event.model_kmer : SequenceUtil.reverseComplement(event.model_kmer);

			prevPos = event.position;
			if (!event.reference_kmer.equals(kmer)) {
				if (match > 0) {
					mis.add(match + "");
					match = 0;
				}
				mis.add(kmer);
			} else
				match++;
		}
		if (match > 0) {
			mis.add(match + "");
			match = 0;
		}

		StringBuffer sb = new StringBuffer();
		for (String m : mis)
			sb.append(m);
		return sb.toString();
	}

	public Cigar calculateCigar(EASeries series) {
		int prevPos = series.position[0];
		ModifyableCigar cigar = new ModifyableCigar();
		cigar.add(series.event_index_start - 1, CigarOperator.H);
		cigar.add(1, CigarOperator.M);
		for (int i = 1; i < series.position.length; i++) {
			int delta = series.position[i] - prevPos;
			switch (delta) {
			case 0:
				cigar.add(1, CigarOperator.I);
				break;
			case 1:
				cigar.add(1, CigarOperator.M);
				break;
			default:
				cigar.add(delta - 1, CigarOperator.D);
				cigar.add(1, CigarOperator.M);
				break;
			}

			// System.out.printf("i=%d, pos=%d, delta=%d, cigar=%s\n", i,
			// series.position[i], delta, cigar);
			prevPos = series.position[i];
		}

		return cigar.toCigar();
	}

	private static class ModifyableCigar {
		LinkedList<ModifyableCigarElement> elements = new LinkedList<EventAlignReader.ModifyableCigarElement>();

		void add(int len, CigarOperator op) {
			if (elements.isEmpty()) {
				elements.add(new ModifyableCigarElement(len, op));
				return;
			}

			ModifyableCigarElement last = elements.getLast();
			if (last.op == op)
				last.len += len;
			else
				elements.add(new ModifyableCigarElement(len, op));
		}

		Cigar toCigar() {
			List<CigarElement> cigarElements = new ArrayList<CigarElement>(elements.size());
			for (ModifyableCigarElement mce : elements)
				cigarElements.add(mce.toCigarElement());

			return new Cigar(cigarElements);
		}

		@Override
		public String toString() {
			return toCigar().toString();
		}
	}

	private static class ModifyableCigarElement {
		int len;
		CigarOperator op;

		public ModifyableCigarElement(int len, CigarOperator op) {
			super();
			this.len = len;
			this.op = op;
		}

		public CigarElement toCigarElement() {
			return new CigarElement(len, op);
		}
	}

	private static class Event {
		String contig;
		int position;
		String reference_kmer;
		int read_index;
		boolean strand;
		int event_index;
		float event_level_mean;
		float event_length;
		String model_kmer;
		float model_mean;
		float model_stdv;
		String model_name;

		public static Event fromString(String line) {
			Event event = new Event();

			try {
				String[] words = line.split("\\s+");
				int index = 0;
				event.contig = words[index++];
				event.position = Integer.parseInt(words[index++]);
				event.reference_kmer = words[index++];
				event.read_index = Integer.parseInt(words[index++]);
				event.strand = "t".equals(words[index++]);
				event.event_index = Integer.parseInt(words[index++]);
				event.event_level_mean = Float.parseFloat(words[index++]);
				event.event_length = Float.parseFloat(words[index++]);
				event.model_kmer = words[index++];
				event.model_mean = Float.parseFloat(words[index++]);
				event.model_stdv = Float.parseFloat(words[index++]);
				event.model_name = words[index++];
			} catch (Exception e) {
				System.err.println("Failed to parse line: " + line);
				e.printStackTrace();
			}

			return event;
		}
	}

	public static void main(String[] args) throws IOException {
		File tsvFile = new File(args[0]);
		File meanSignalFile = new File(args[1]);

		EventAlignReader reader = new EventAlignReader(tsvFile);
		if (reader.firstEvent == null)
			return;
		List<Event> nextRead = null;
		SAMFileHeader header = new SAMFileHeader();
		header.addSequence(new SAMSequenceRecord(reader.firstEvent.contig));
		SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, false,
				new File(tsvFile.getAbsolutePath() + ".bam"));

		int refLen = 10 * 1024 * 1024;
		byte[] ref = new byte[refLen];
		Arrays.fill(ref, (byte) 'N');
		float[] signalT = new float[refLen];
		short[] depthT = new short[refLen];
		float[] signalC = new float[refLen];
		short[] depthC = new short[refLen];
		int maxRefPos = 0;
		long counter = 0;
		while ((nextRead = reader.nextRead()) != null) {
			counter++;
			EASeries series = reader.transpose(nextRead);
			for (int i = 0; i < nextRead.size(); i++) {
				Event e = nextRead.get(i);
				int pos = e.position;
				float[] signal = e.strand ? signalT : signalC;
				short[] depth = e.strand ? depthT : depthC;
				signal[pos] += e.event_level_mean;
				depth[pos]++;
				if (maxRefPos < pos)
					maxRefPos = pos;
			}

			Cigar cigar = reader.calculateCigar(series);
			byte[] bases = reader.calculateBases(nextRead);
			byte[] scores = reader.calculateScores(nextRead);
			String ed = reader.calculateED(nextRead);
			String md = reader.calculateMD(nextRead);
			// System.out.printf("%s\t%s\t%s\tED:Z:%s\tMD:Z:%s\n",
			// cigar.toString(), StringUtil.bytesToString(bases),
			// SAMUtils.phredToFastq(scores), ed, md);

			Event event = nextRead.get(0);
			SAMRecord record = new SAMRecord(header);
			record.setReadName(event.read_index + "");
			record.setReferenceName(event.contig);
			record.setReadNegativeStrandFlag(!event.strand);
			record.setAlignmentStart(5 + event.position);
			record.setCigar(cigar);
			record.setReadBases(bases);
			// record.setReadBases(SAMRecord.NULL_SEQUENCE);
			record.setBaseQualities(scores);
			record.setAttribute("ED", ed);
			// record.setAttribute("MK", series.model_kmer);
			record.setAttribute("EL", series.event_length);
			record.setAttribute("EM", series.event_level_mean);
			// record.setAttribute("ER", md);

			writer.addAlignment(record);

			for (AlignmentBlock b : record.getAlignmentBlocks()) {
				for (int i = 0; i < b.getLength(); i++) {
					ref[b.getReferenceStart() + i - 1] = bases[i + b.getReadStart() - 1];
				}
			}
		}
		writer.close();

		if (counter > 0) {
			PrintStream ps = new PrintStream(meanSignalFile);
			for (int i = 0; i < maxRefPos; i++) {
				ps.printf("%.2f\t%.2f\n", signalT[i] / depthT[i], signalC[i] / depthC[i]);
			}
			ps.close();
		}
	}

	private static SAMRecord createRecord(SAMFileHeader header, String[] data) {
		SAMRecord record = createNewRecord(data[0], header);
		EASeries s = getSignal(data);
		record.setAttribute("EL", s.event_length);
		record.setAttribute("EM", s.event_level_mean);
		record.setAttribute("MK", s.model_kmer);

		List<Byte> bases = new ArrayList<Byte>();
		{
			int prev = 0;
			List<String> mis = new ArrayList<String>();
			// bases.add(kmers[s.model_kmer[0]].getBytes()[0]);
			// bases.add(kmers[s.model_kmer[0]].getBytes()[1]);
			// bases.add(kmers[s.model_kmer[0]].getBytes()[2]);
			// bases.add(kmers[s.model_kmer[0]].getBytes()[3]);
			prev = s.model_kmer[0] >> 2;
			int i = 0;
			for (int c : s.model_kmer) {
				byte base;
				if ((c >> 2) == (prev & 255)) {
					mis.add(EventAlign.kmers[c]);
					base = ((byte) 'N');
					base = EventAlign.kmers[c].getBytes()[4];
					bases.add(base);
				} else if ((c >> 4) == (prev & 63)) {
					mis.add(EventAlign.kmers[c]);
					base = ((byte) 'N');
					base = EventAlign.kmers[c].getBytes()[3];
					bases.add(base);
					base = EventAlign.kmers[c].getBytes()[4];
					bases.add(base);
				} else if ((c >> 6) == (prev & 15)) {
					mis.add(EventAlign.kmers[c]);
					base = ((byte) 'N');
					base = EventAlign.kmers[c].getBytes()[2];
					bases.add(base);
					base = EventAlign.kmers[c].getBytes()[3];
					bases.add(base);
					base = EventAlign.kmers[c].getBytes()[4];
					bases.add(base);
				} else {
					base = ((byte) 'N');
					bases.add(base);
				}

				// System.out.println(data[i]);
				// System.out.printf("i=%d: c=%d, kmer=%s, base=%c\n", i, c,
				// kmers[c], (char) base);
				// if (c != prev)

				prev = c;
				i++;
			}
		}

		byte[] b = new byte[bases.size()];
		// bases.set(0, (byte) '!');
		for (int i = 0; i < b.length; i++) {
			b[i] = bases.get(i);
		}

		record.setReadBases(b);
		record.setBaseQualities(s.scores);
		record.setReadName("" + s.read_index);

		List<CigarElement> l = new ArrayList<CigarElement>();
		if (s.event_index_start > 1)
			l.add(new CigarElement(s.event_index_start - 1, CigarOperator.H));
		int pos = s.position[0];
		int match = 1, ins = 0, del = 0;

		for (int i = 1; i < s.position.length; i++) {
			int len = s.position[i] - pos;
			if (len == 1) {
				if (ins != 0)
					l.add(new CigarElement(ins, CigarOperator.I));
				if (del != 0)
					l.add(new CigarElement(del, CigarOperator.D));
				match++;
				ins = 0;
				del = 0;
			} else if (len == 0) {
				if (match != 0)
					l.add(new CigarElement(match, CigarOperator.M));
				if (del != 0)
					l.add(new CigarElement(del, CigarOperator.D));
				ins++;
				match = 0;
				del = 0;
			} else {
				if (match != 0)
					l.add(new CigarElement(match, CigarOperator.M));
				if (ins != 0)
					l.add(new CigarElement(ins, CigarOperator.I));
				del++;
				match = 0;
				ins = 0;
			}
			pos = s.position[i];
		}
		if (match != 0)
			l.add(new CigarElement(match, CigarOperator.M));
		if (ins != 0)
			l.add(new CigarElement(ins, CigarOperator.I));
		if (del != 0)
			l.add(new CigarElement(del, CigarOperator.D));

		if (new Cigar(l).getReadLength() < s.position.length) {
			if (l.get(l.size() - 1).getOperator() == CigarOperator.M) {
				int len = l.remove(l.size() - 1).getLength() + s.position.length - new Cigar(l).getReadLength();
				l.add(new CigarElement(len, CigarOperator.M));
			} else
				l.add(new CigarElement(s.position.length - new Cigar(l).getReadLength(), CigarOperator.M));

		}
		record.setCigar(new Cigar(l));

		return record;
	}

	private static class EASeries {
		float[] event_level_mean;
		float[] event_length;
		int[] model_kmer;
		int[] position;
		int event_index_start;
		int read_index;
		boolean template;
		byte[] scores;

		public EASeries(int size) {
			event_length = new float[size];
			event_level_mean = new float[size];
			model_kmer = new int[size];
			position = new int[size];
			scores = new byte[size];
		}
	}

	private static EASeries getSignal(String[] data) {
		EASeries s = new EASeries(data.length);
		int i = 0;
		for (String line : data) {
			String[] words = line.split("\\s+");
			if (i == 0) {
				s.event_index_start = Integer.parseInt(words[EA.event_index.ordinal()]);
				s.read_index = Integer.parseInt(words[EA.read_index.ordinal()]);
				s.template = "t".equals(words[EA.read_index.ordinal()]);
			}
			s.event_length[i] = Float.parseFloat(words[EA.event_length.ordinal()]);
			s.event_level_mean[i] = Float.parseFloat(words[EA.event_level_mean.ordinal()]);
			String kmer = words[EA.model_kmer.ordinal()];
			if (words[EA.strand.ordinal()].equals("c")) {
				kmer = SequenceUtil.reverseComplement(kmer);
				System.out.println(kmer);
			}
			int n = 0;

			for (byte b : kmer.getBytes()) {
				n <<= 2;
				n = n | EventAlign.m[b];
			}
			s.model_kmer[i] = EventAlign.codes.get(kmer);
			if (!EventAlign.kmers[n].equals(kmer)) {
				System.err.printf("kmer=%s, n=%d, kmersn=%s\n", kmer, n, EventAlign.kmers[n]);
				throw new RuntimeException();
			}

			s.position[i] = Integer.parseInt(words[EA.position.ordinal()]);
			float elm = Float.parseFloat(words[EA.event_level_mean.ordinal()]);
			float mm = Float.parseFloat(words[EA.model_mean.ordinal()]);
			float stdv = Float.parseFloat(words[EA.model_stdv.ordinal()]);
			s.scores[i] = score(elm, mm, stdv);
			// System.out.printf("elm=%.2f, mm=%.2f, stdv=%.2f, s=%c, s=%d\n",
			// elm, mm, stdv, (char) s.scores[i], s.scores[i]);
			i++;
		}
		return s;
	}

	private static SAMRecord createNewRecord(String ea, SAMFileHeader header) {
		String[] words = ea.split("\\s+");
		SAMRecord record = new SAMRecord(header);
		record.setReferenceName(words[EA.contig.ordinal()]);
		record.setAlignmentStart(Integer.parseInt(words[EA.position.ordinal()]));
		record.setReadNegativeStrandFlag(!"t".equals(words[EA.strand.ordinal()]));
		record.setAttribute("EI", words[EA.event_index.ordinal()]);
		return record;
	}

	private static byte score(float x, float mean, float stdv) {
		float f = normalDistr(x, mean, stdv);
		return (byte) (f * 40);
	}

	private static float normalDistr(float x, float mean, float stdv) {
		return (float) (1f / Math.sqrt(2f * Math.PI) / stdv * Math.exp(-Math.pow(x - mean, 2) / 2f / stdv / stdv));
	}
}
