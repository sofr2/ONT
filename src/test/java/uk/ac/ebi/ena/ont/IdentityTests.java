package uk.ac.ebi.ena.ont;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMRecord;

import org.junit.Test;

public class IdentityTests {

	@Test
	public void test_2M() {
		Identity id = new Identity("AAAC".getBytes());
		SAMFileHeader header = new SAMFileHeader();
		SAMRecord r = new SAMRecord(header);
		r.setAlignmentStart(1);
		r.setReadBases("AA".getBytes());
		r.setCigarString("2M");
		int identical = id.addRead(r);
		assertThat(identical, is(2));
		assertThat(id.getAlignedBaseCount(), is(2L));
		assertThat(id.getIdenticalBaseCount(), is(2L));
		assertThat(id.getNonIdenticalBaseCount(), is(0L));
		assertThat(id.getTotalBaseCount(), is(2L));
		assertThat(id.getUnalignedBaseCount(), is(0L));
		assertThat(id.getUnmappedBaseCount(), is(0L));
	}

	@Test
	public void test_2M1I() {
		Identity id = new Identity("AAAC".getBytes());
		SAMFileHeader header = new SAMFileHeader();
		SAMRecord r = new SAMRecord(header);
		r.setAlignmentStart(1);
		r.setReadBases("AAT".getBytes());
		r.setCigarString("2M1I");
		int identical = id.addRead(r);
		assertThat(identical, is(2));
		assertThat(id.getAlignedBaseCount(), is(2L));
		assertThat(id.getIdenticalBaseCount(), is(2L));
		assertThat(id.getNonIdenticalBaseCount(), is(0L));
		assertThat(id.getTotalBaseCount(), is(3L));
		assertThat(id.getUnalignedBaseCount(), is(1L));
		assertThat(id.getUnmappedBaseCount(), is(0L));
	}

	@Test
	public void test_2M1D() {
		Identity id = new Identity("AAAC".getBytes());
		SAMFileHeader header = new SAMFileHeader();
		SAMRecord r = new SAMRecord(header);
		r.setAlignmentStart(1);
		r.setReadBases("AA".getBytes());
		r.setCigarString("2M1D");
		int identical = id.addRead(r);
		assertThat(identical, is(2));
		assertThat(id.getAlignedBaseCount(), is(2L));
		assertThat(id.getIdenticalBaseCount(), is(2L));
		assertThat(id.getNonIdenticalBaseCount(), is(0L));
		assertThat(id.getTotalBaseCount(), is(2L));
		assertThat(id.getUnalignedBaseCount(), is(0L));
		assertThat(id.getUnmappedBaseCount(), is(0L));
	}

	@Test
	public void test_1match_1mismatch() {
		Identity id = new Identity("AAAC".getBytes());
		SAMFileHeader header = new SAMFileHeader();
		SAMRecord r = new SAMRecord(header);
		r.setAlignmentStart(1);
		r.setReadBases("AT".getBytes());
		r.setCigarString("2M");
		int identical = id.addRead(r);
		assertThat(identical, is(1));
		assertThat(id.getAlignedBaseCount(), is(2L));
		assertThat(id.getIdenticalBaseCount(), is(1L));
		assertThat(id.getNonIdenticalBaseCount(), is(1L));
		assertThat(id.getTotalBaseCount(), is(2L));
		assertThat(id.getUnalignedBaseCount(), is(0L));
		assertThat(id.getUnmappedBaseCount(), is(0L));
	}

	@Test
	public void test_1mapped_1unmapped() {
		Identity id = new Identity("AAAC".getBytes());
		SAMFileHeader header = new SAMFileHeader();
		SAMRecord r1 = new SAMRecord(header);
		r1.setAlignmentStart(1);
		r1.setReadBases("AA".getBytes());
		r1.setCigarString("2M");
		int identical = id.addRead(r1);
		assertThat(identical, is(2));

		SAMRecord r2 = new SAMRecord(header);
		r2.setReadBases("NNN".getBytes());
		r2.setReadUnmappedFlag(true);
		identical = id.addRead(r2);
		assertThat(identical, is(0));

		assertThat(id.getAlignedBaseCount(), is(2L));
		assertThat(id.getIdenticalBaseCount(), is(2L));
		assertThat(id.getNonIdenticalBaseCount(), is(0L));
		assertThat(id.getTotalBaseCount(), is(5L));
		assertThat(id.getUnalignedBaseCount(), is(0L));
		assertThat(id.getUnmappedBaseCount(), is(3L));
	}

}
