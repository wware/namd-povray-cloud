/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef MSMMAP_H
#define MSMMAP_H

#define DEBUG_MSM
#undef DEBUG_MSM

#define DEBUG_MSM_VERBOSE
#undef DEBUG_MSM_VERBOSE

#define DEBUG_MSM_GRID
#undef DEBUG_MSM_GRID

// assert macro
#undef ASSERT
#ifdef DEBUG_MSM
#define ASSERT(expr) \
  do { \
    if ( !(expr) ) { \
      char msg[100]; \
      snprintf(msg, sizeof(msg), "ASSERT: \"%s\" " \
          "(%s, %d)\n", #expr, __FILE__, __LINE__); \
      NAMD_die(msg); \
    } \
  } while (0)
#else
#define ASSERT(expr)
#endif 


namespace msm {

  ///////////////////////////////////////////////////////////////////////////
  //
  // Resizable Array class
  //
  ///////////////////////////////////////////////////////////////////////////

  template <class T> class Array;

  template <class T>
  void swap(Array<T>& s, Array<T>& t);

  template <class T>
  class Array {
    public:
      Array() : abuffer(0), alen(0), amax(0) { }
      Array(int n) : abuffer(0), alen(0), amax(0) { resize(n); }
      Array(const Array& a) : abuffer(0), alen(0), amax(0) { copy(a); }
      ~Array() { resize(0); }
      Array& operator=(const Array& a) {
        if (this != &a) copy(a);  // don't allow self-assignment
        return(*this);
      }
      int len() const { return alen; }
      int max() const { return amax; }
      const T& operator[](int i) const {
#ifdef DEBUG_MSM
        return elem(i);
#else
        return abuffer[i];
#endif
      }
      const T& elem(int i) const {
        if (i < 0 || i >= alen) {
          char msg[100];
          snprintf(msg, sizeof(msg), "Array index:  alen=%d, i=%d\n", alen, i);
          NAMD_die(msg);
        }
        return abuffer[i];
      }
      T& operator[](int i) {
#ifdef DEBUG_MSM
        return elem(i);
#else
        return abuffer[i];
#endif
      }
      T& elem(int i) {
        if (i < 0 || i >= alen) {
          char msg[100];
          snprintf(msg, sizeof(msg), "Array index:  alen=%d, i=%d\n", alen, i);
          NAMD_die(msg);
        }
        return abuffer[i];
      }
      void append(const T& t) {
        if (alen==amax) setmax(2*amax+1);
        abuffer[alen++] = t;
      }
      void resize(int n) {
        if (n > amax) setmax(n);
        alen = n;
      }
      void setmax(int m);
      const T *buffer() const { return abuffer; }
      T *buffer() { return abuffer; }
      const T *buffer(int& n) const { n = alen; return abuffer; }
      T *buffer(int& n) { n = alen; return abuffer; }
      friend void swap<T>(Array&, Array&);
#ifdef DEBUG_MSM
      void print(const char *s=0) const {
        if (s) printf("PRINTING DATA FOR ARRAY \"%s\":\n", s);
        printf("abuffer=%p\n  alen=%d  amax=%d\n",
            (void *) abuffer, alen, amax);
      }
#endif
    protected:
      T *abuffer;
      int alen, amax;
      void copy(const Array& a);
  };

  template <class T>
  void Array<T>::setmax(int m) {
    if (m == amax) return;
    else if (m > 0) {
      T *newbuffer = new T[m];
      if ( ! newbuffer) {
        char msg[100];
        snprintf(msg, sizeof(msg), "Can't allocate %lu KB\n",
            (unsigned long)(m * sizeof(T) / 1024));
        NAMD_die(msg);
      }
      if (alen > m) alen = m;  // new buffer is shorter than old buffer
      for (int i = 0;  i < alen;  i++) {
        newbuffer[i] = abuffer[i];
      }
      delete[] abuffer;
      abuffer = newbuffer;
      amax = m;
    }
    else {  // consider m == 0
      delete[] abuffer;
      abuffer = 0;
      alen = 0;
      amax = 0;
    }
  }

  template <class T>
  void Array<T>::copy(const Array<T>& a) {
    setmax(a.amax);
    alen = a.alen;
    for (int i = 0;  i < alen;  i++) {
      abuffer[i] = a.abuffer[i];
    }
  }

  // swap arrays without duplicating memory buffer
  template <class T>
  void swap(Array<T>& s, Array<T>& t) {
    T *tmpbuffer = s.abuffer;  s.abuffer = t.abuffer;  t.abuffer = tmpbuffer;
    tmpbuffer = 0;
    int tmpn = s.alen;  s.alen = t.alen;  t.alen = tmpn;
    tmpn = s.amax;  s.amax = t.amax;  t.amax = tmpn;
    tmpn = s.astate;  s.astate = t.astate;  t.astate = tmpn;
  }


  ///////////////////////////////////////////////////////////////////////////
  //
  // Grid is 3D lattice of grid points with user-definable index ranges.
  //
  ///////////////////////////////////////////////////////////////////////////

  // 3-integer vector, used for indexing from a 3D grid
  struct Ivec {
    Ivec(int n=0) : i(n), j(n), k(n) { }
    Ivec(int ni, int nj, int nk) : i(ni), j(nj), k(nk) { }
    int operator==(const Ivec& n) { return(i==n.i && j==n.j && k==n.k); }
    int i, j, k;
  };

  // index range for 3D lattice of grid points
  class IndexRange {
    public:
      IndexRange() : nlower(0), nextent(0) { }
      void set(int pia, int pni, int pja, int pnj, int pka, int pnk) {
        ASSERT(pni >= 0 && pnj >= 0 && pnk >= 0);
        nlower = Ivec(pia, pja, pka);
        nextent = Ivec(pni, pnj, pnk);
      }
      void setbounds(int pia, int pib, int pja, int pjb, int pka, int pkb) {
        set(pia, pib-pia+1, pja, pjb-pja+1, pka, pkb-pka+1);
      }
      int ia() const { return nlower.i; }
      int ib() const { return nlower.i + nextent.i - 1; }
      int ja() const { return nlower.j; }
      int jb() const { return nlower.j + nextent.j - 1; }
      int ka() const { return nlower.k; }
      int kb() const { return nlower.k + nextent.k - 1; }
      int ni() const { return nextent.i; }
      int nj() const { return nextent.j; }
      int nk() const { return nextent.k; }
      int nn() const { return nextent.i * nextent.j * nextent.k; }
      Ivec lower() const { return nlower; }
      Ivec extent() const { return nextent; }
      int operator<=(const IndexRange& n) {
        // true if this IndexRange fits inside n
        return ( ia() >= n.ia() && ib() <= n.ib() &&
                 ja() >= n.ja() && jb() <= n.jb() &&
                 ka() >= n.ka() && kb() <= n.kb() );
      }
    protected:
      Ivec nlower;   // index for lowest corner of rectangular lattice
      Ivec nextent;  // extent of lattice along each dimension
  };

  // storage and indexing for 3D lattice of grid points
  template <class T>
  class Grid : public IndexRange {
    public:
      Grid() { }
      void init(const IndexRange& n) {
        nlower = n.lower();
        nextent = n.extent();
        gdata.resize(nn());
      }
      void set(int pia, int pni, int pja, int pnj, int pka, int pnk) {
        IndexRange::set(pia, pni, pja, pnj, pka, pnk);
        gdata.resize(nn());
      }
      void setbounds(int pia, int pib, int pja, int pjb, int pka, int pkb) {
        IndexRange::setbounds(pia, pib, pja, pjb, pka, pkb);
        gdata.resize(nn());
      }
      void resize(int n) { // reserve space but don't set grid indexing
        gdata.resize(n);
      }
      const T& operator()(int i, int j, int k) const {
#ifdef DEBUG_MSM
        return elem(i,j,k);
#else
        return gdata[flatindex(i,j,k)];
#endif
      }
      const T& operator()(const Ivec& n) const {
        return this->operator()(n.i, n.j, n.k);
      }
      const T& elem(int i, int j, int k) const {
        if (i<ia() || i>ib() || j<ja() || j>jb() || k<ka() || k>kb()) {
          char msg[200];
          snprintf(msg, sizeof(msg), "Grid indexing:\n"
              "ia=%d, ib=%d, i=%d\n"
              "ja=%d, jb=%d, j=%d\n"
              "ka=%d, kb=%d, k=%d\n",
              ia(), ib(), i, ja(), jb(), j, ka(), kb(), k);
          NAMD_die(msg);
        }
        return gdata[flatindex(i,j,k)];
      }
      T& operator()(int i, int j, int k) {
#ifdef DEBUG_MSM
        return elem(i,j,k);
#else
        return gdata[flatindex(i,j,k)];
#endif
      }
      T& operator()(const Ivec& n) {
        return this->operator()(n.i, n.j, n.k);
      }
      T& elem(int i, int j, int k) {
        if (i<ia() || i>ib() || j<ja() || j>jb() || k<ka() || k>kb()) {
          char msg[200];
          snprintf(msg, sizeof(msg), "Grid indexing:\n"
              "ia=%d, ib=%d, i=%d\n"
              "ja=%d, jb=%d, j=%d\n"
              "ka=%d, kb=%d, k=%d\n",
              ia(), ib(), i, ja(), jb(), j, ka(), kb(), k);
          NAMD_die(msg);
        }
        return gdata[flatindex(i,j,k)];
      }
      int flatindex(int i, int j, int k) const {
        return ((k-ka())*nj() + (j-ja()))*ni() + (i-ia());
      }
      const Array<T>& data() const { return gdata; }
      Array<T>& data() { return gdata; }

      // use to zero out grid
      void reset(const T& t) {
        for (int n = 0;  n < gdata.len();  n++) { gdata[n] = t; }
      }

      // use to modify the indexing by changing lower corner
      void updateLower(const Ivec& n) { nlower = n; }

      // accumulate another grid into this grid
      // the grid to be added must fit within this grid's index range
      Grid<T>& operator+=(const Grid<T>& g) {
        ASSERT(IndexRange(g) <= IndexRange(*this));
        int gia = g.ia();
        int gib = g.ib();
        int gja = g.ja();
        int gjb = g.jb();
        int gka = g.ka();
        int gkb = g.kb();
        for (int k = gka;  k <= gkb;  k++) {
          for (int j = gja;  j <= gjb;  j++) {
            for (int i = gia;  i <= gib;  i++) {
              (*this)(i,j,k) += g(i,j,k);
            }
          }
        }
        return(*this);
      }

      // extract a subgrid from this grid
      // subgrid must fit within this grid's index range
      void extract(Grid<T>& g) {
        ASSERT(IndexRange(g) <= IndexRange(*this));
        int gia = g.ia();
        int gib = g.ib();
        int gja = g.ja();
        int gjb = g.jb();
        int gka = g.ka();
        int gkb = g.kb();
        for (int k = gka;  k <= gkb;  k++) {
          for (int j = gja;  j <= gjb;  j++) {
            for (int i = gia;  i <= gib;  i++) {
              g(i,j,k) = (*this)(i,j,k);
            }
          }
        }
      }

    private:
      Array<T> gdata;
  };


  ///////////////////////////////////////////////////////////////////////////
  //
  // Map object 
  //
  ///////////////////////////////////////////////////////////////////////////

  // index a block from the MSM grid hierarchy
  struct BlockIndex {
    int level;
    Ivec n;
    BlockIndex() : level(0), n(0) { }
    BlockIndex(int ll, const Ivec& nn) : level(ll), n(nn) { }
  };

  // sending part of an extended grid calculation to another block
  struct BlockSend {
    BlockIndex nblock;       // relative block index
    IndexRange nrange;       // relative grid index range
    BlockIndex nblock_wrap;  // true block index
    IndexRange nrange_wrap;  // true grid index range
    void reset() {
      nblock = BlockIndex();
      nrange = IndexRange();
      nblock_wrap = BlockIndex();
      nrange_wrap = IndexRange();
    } // reset
  };

  struct PatchSend {
    IndexRange nrange;         // true grid index range from my block
    IndexRange nrange_unwrap;  // relative grid index range for patch
    int patchID;               // patch ID
    void reset() {
      nrange = IndexRange();
      nrange_unwrap = IndexRange();
      patchID = -1;
    } // reset
  };

  // one PatchDiagram for each patch
  // maintain a Grid of PatchDiagram, indexed by patch ID
  struct PatchDiagram {
    IndexRange nrange;       // shows subset of MSM h-grid covering this patch
    Array<BlockSend> send;   // array of blocks to which this patch sends
    int numRecvs;            // number of blocks from which this patch receives
    void reset() {
      nrange = IndexRange();
      send.resize(0);
      numRecvs = 0;
    } // reset
  };

  // one BlockDiagram for each block of each level of each MSM grid
  // maintain a Grid of BlockDiagram for each level
  struct BlockDiagram {
    IndexRange nrange;            // subset of MSM grid for this block
    IndexRange nrangeCutoff;      // expanded subgrid for cutoff calculation
    IndexRange nrangeRestricted;  // (level+1) subgrid for restriction
    IndexRange nrangeProlongated; // (level-1) subgrid for prolongation
    Array<BlockSend> sendUp;      // send up charge to blocks on (level+1)
    Array<BlockSend> sendAcross;  // send across potential to blocks on (level)
    Array<BlockSend> sendDown;    // send down potential to blocks on (level-1)
    Array<PatchSend> sendPatch;   // send my (level=0) potential block to patch
    int numRecvsCharge;           // number of expected receives of charge
    int numRecvsPotential;        // number of expected receives of potential
    void reset() {
      nrange = IndexRange();
      nrangeCutoff = IndexRange();
      nrangeRestricted = IndexRange();
      nrangeProlongated = IndexRange();
      sendUp.resize(0);
      sendAcross.resize(0);
      sendDown.resize(0);
      sendPatch.resize(0);
      numRecvsCharge = 0;
      numRecvsPotential = 0;
    } // reset
  };


  struct Map {
    Array<IndexRange> gridrange;  // dimensions for each MSM grid level

    Array<Grid<BigReal> > gc;     // grid constant weights for each level

    Array<PatchDiagram> patchList;
    Array<Grid<BlockDiagram> > blockLevel;

    int ispx, ispy, ispz;         // is periodic in x, y, z?

    Array<int> bsx, bsy, bsz;     // block size in x, y, z for each level

    // clip index to grid level, using periodicity flags
    Ivec clipIndexToLevel(const Ivec& n, int level) const {
      ASSERT(level >= 0 && level < gridrange.len());
      Ivec pn(n);
      if ( ! ispx) {
        int a = gridrange[level].ia();
        int b = gridrange[level].ib();
        if (pn.i < a) pn.i = a;
        if (pn.i > b) pn.i = b;
      }
      if ( ! ispy) {
        int a = gridrange[level].ja();
        int b = gridrange[level].jb();
        if (pn.j < a) pn.j = a;
        if (pn.j > b) pn.j = b;
      }
      if ( ! ispz) {
        int a = gridrange[level].ka();
        int b = gridrange[level].kb();
        if (pn.k < a) pn.k = a;
        if (pn.k > b) pn.k = b;
      }
      return pn;
    }

    // determine relative (unwrapped) block index for the given grid index
    BlockIndex blockOfGridIndex(const Ivec& n, int level) const {
      ASSERT(level >= 0 && level < gridrange.len());
      BlockIndex bn;
      // we want floor((i - ia) / bsx), etc.
      // modify case i < ia to avoid integer division of negative numbers
      int d = n.i - gridrange[level].ia();
      bn.n.i = (d >= 0 ? d / bsx[level] : -((-d+bsx[level]-1) / bsx[level]));
      d = n.j - gridrange[level].ja();
      bn.n.j = (d >= 0 ? d / bsy[level] : -((-d+bsy[level]-1) / bsy[level]));
      d = n.k - gridrange[level].ka();
      bn.n.k = (d >= 0 ? d / bsz[level] : -((-d+bsz[level]-1) / bsz[level]));
      bn.level = level;
      return bn;
    }

    // find the natural index range of the given relative block number
    IndexRange indexRangeOfBlock(const BlockIndex& nb) const {
      ASSERT(nb.level >= 0 && nb.level < gridrange.len());
      IndexRange nr;
      int ia = nb.n.i * bsx[nb.level] + gridrange[nb.level].ia();
      int ja = nb.n.j * bsy[nb.level] + gridrange[nb.level].ja();
      int ka = nb.n.k * bsz[nb.level] + gridrange[nb.level].ka();
      nr.set(ia, bsx[nb.level], ja, bsy[nb.level], ka, bsz[nb.level]);
      return nr;
    }

    // clip the natural block index range to not exceed the given index range
    IndexRange clipBlockToIndexRange(const BlockIndex& nb,
        const IndexRange& nrange) const {
      IndexRange nr = indexRangeOfBlock(nb);
      int nia = nrange.ia();
      int nib = nrange.ib();
      int nja = nrange.ja();
      int njb = nrange.jb();
      int nka = nrange.ka();
      int nkb = nrange.kb();
      int ia = nr.ia();
      if (ia < nia) ia = nia;
      int ib = nr.ib();
      if (ib > nib) ib = nib;
      int ja = nr.ja();
      if (ja < nja) ja = nja;
      int jb = nr.jb();
      if (jb > njb) jb = njb;
      int ka = nr.ka();
      if (ka < nka) ka = nka;
      int kb = nr.kb();
      if (kb > nkb) kb = nkb;
      nr.setbounds(ia, ib, ja, jb, ka, kb);
      return nr;
    }

    // set the nblock_wrap and nrange_wrap fields based on periodicity
    void wrapBlockSend(BlockSend& bs) const {
      BlockIndex nb = bs.nblock;
      IndexRange nr = bs.nrange;
      int level = bs.nblock.level;
      ASSERT(level >= 0 && level < blockLevel.len());
      int ni = blockLevel[level].ni();
      int nj = blockLevel[level].nj();
      int nk = blockLevel[level].nk();
      int di=0, dj=0, dk=0;
      if (ispx) {
        while (nb.n.i < 0) {
          nb.n.i += ni;
          di += ni * bsx[level];
        }
        while (nb.n.i >= ni) {
          nb.n.i -= ni;
          di -= ni * bsx[level];
        }
      }
      if (ispy) {
        while (nb.n.j < 0) {
          nb.n.j += nj;
          dj += nj * bsy[level];
        }
        while (nb.n.j >= nj) {
          nb.n.j -= nj;
          dj -= nj * bsy[level];
        }
      }
      if (ispz) {
        while (nb.n.k < 0) {
          nb.n.k += nk;
          dk += nk * bsz[level];
        }
        while (nb.n.k >= nk) {
          nb.n.k -= nk;
          dk -= nk * bsz[level];
        }
      }
      int ia = nr.ia();
      int ib = nr.ib();
      int ja = nr.ja();
      int jb = nr.jb();
      int ka = nr.ka();
      int kb = nr.kb();
      nr.setbounds(ia + di, ib + di, ja + dj, jb + dj, ka + dk, kb + dk);
      bs.nblock_wrap = nb;
      bs.nrange_wrap = nr;
    }
  }; // Map


  struct AtomCoord {
    Position position;
    float charge;
    int id;
  };

  typedef Array<AtomCoord> AtomCoordArray;
  typedef Array<Force> ForceArray;

  struct PatchData;
  typedef Array<PatchData> PatchDataArray;
  typedef Array<PatchData *> PatchPtrArray;

  struct BlockData;
  typedef Array<Grid<BlockData> > BlockDataGrids;

} // namespace msm

#endif // MSMMAP_H
