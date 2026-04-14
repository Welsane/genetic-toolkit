/**
 * app.js — Genetic Analyzer Toolkit Frontend
 * All bioinformatics logic runs client-side in JS, mirroring the C backend.
 * To connect to a real C backend, replace the compute functions with fetch()
 * calls to your API endpoint (see API_BASE constant below).
 */

'use strict';

/* ── API Config (swap URL to connect real C backend) ───────────────── */
const USE_API   = false;          // set true + define URL to use C backend
const API_BASE  = 'http://localhost:8080/api';

/* ═══════════════════════════════════════════════════════════════════════
   CODON TABLE  (standard genetic code, DNA codons)
═══════════════════════════════════════════════════════════════════════ */
const CODON_TABLE = {
  TTT:{name:'Phenylalanine',abbr:'Phe',letter:'F',type:'nonpolar'},
  TTC:{name:'Phenylalanine',abbr:'Phe',letter:'F',type:'nonpolar'},
  TTA:{name:'Leucine',abbr:'Leu',letter:'L',type:'nonpolar'},
  TTG:{name:'Leucine',abbr:'Leu',letter:'L',type:'nonpolar'},
  CTT:{name:'Leucine',abbr:'Leu',letter:'L',type:'nonpolar'},
  CTC:{name:'Leucine',abbr:'Leu',letter:'L',type:'nonpolar'},
  CTA:{name:'Leucine',abbr:'Leu',letter:'L',type:'nonpolar'},
  CTG:{name:'Leucine',abbr:'Leu',letter:'L',type:'nonpolar'},
  ATT:{name:'Isoleucine',abbr:'Ile',letter:'I',type:'nonpolar'},
  ATC:{name:'Isoleucine',abbr:'Ile',letter:'I',type:'nonpolar'},
  ATA:{name:'Isoleucine',abbr:'Ile',letter:'I',type:'nonpolar'},
  ATG:{name:'Methionine',abbr:'Met',letter:'M',type:'start'},
  GTT:{name:'Valine',abbr:'Val',letter:'V',type:'nonpolar'},
  GTC:{name:'Valine',abbr:'Val',letter:'V',type:'nonpolar'},
  GTA:{name:'Valine',abbr:'Val',letter:'V',type:'nonpolar'},
  GTG:{name:'Valine',abbr:'Val',letter:'V',type:'nonpolar'},
  TCT:{name:'Serine',abbr:'Ser',letter:'S',type:'polar'},
  TCC:{name:'Serine',abbr:'Ser',letter:'S',type:'polar'},
  TCA:{name:'Serine',abbr:'Ser',letter:'S',type:'polar'},
  TCG:{name:'Serine',abbr:'Ser',letter:'S',type:'polar'},
  AGT:{name:'Serine',abbr:'Ser',letter:'S',type:'polar'},
  AGC:{name:'Serine',abbr:'Ser',letter:'S',type:'polar'},
  CCT:{name:'Proline',abbr:'Pro',letter:'P',type:'special'},
  CCC:{name:'Proline',abbr:'Pro',letter:'P',type:'special'},
  CCA:{name:'Proline',abbr:'Pro',letter:'P',type:'special'},
  CCG:{name:'Proline',abbr:'Pro',letter:'P',type:'special'},
  ACT:{name:'Threonine',abbr:'Thr',letter:'T',type:'polar'},
  ACC:{name:'Threonine',abbr:'Thr',letter:'T',type:'polar'},
  ACA:{name:'Threonine',abbr:'Thr',letter:'T',type:'polar'},
  ACG:{name:'Threonine',abbr:'Thr',letter:'T',type:'polar'},
  GCT:{name:'Alanine',abbr:'Ala',letter:'A',type:'nonpolar'},
  GCC:{name:'Alanine',abbr:'Ala',letter:'A',type:'nonpolar'},
  GCA:{name:'Alanine',abbr:'Ala',letter:'A',type:'nonpolar'},
  GCG:{name:'Alanine',abbr:'Ala',letter:'A',type:'nonpolar'},
  TAT:{name:'Tyrosine',abbr:'Tyr',letter:'Y',type:'polar'},
  TAC:{name:'Tyrosine',abbr:'Tyr',letter:'Y',type:'polar'},
  TAA:{name:'STOP',abbr:'Stop',letter:'*',type:'stop'},
  TAG:{name:'STOP',abbr:'Stop',letter:'*',type:'stop'},
  TGA:{name:'STOP',abbr:'Stop',letter:'*',type:'stop'},
  CAT:{name:'Histidine',abbr:'His',letter:'H',type:'charged'},
  CAC:{name:'Histidine',abbr:'His',letter:'H',type:'charged'},
  CAA:{name:'Glutamine',abbr:'Gln',letter:'Q',type:'polar'},
  CAG:{name:'Glutamine',abbr:'Gln',letter:'Q',type:'polar'},
  AAT:{name:'Asparagine',abbr:'Asn',letter:'N',type:'polar'},
  AAC:{name:'Asparagine',abbr:'Asn',letter:'N',type:'polar'},
  AAA:{name:'Lysine',abbr:'Lys',letter:'K',type:'charged'},
  AAG:{name:'Lysine',abbr:'Lys',letter:'K',type:'charged'},
  GAT:{name:'Aspartate',abbr:'Asp',letter:'D',type:'charged'},
  GAC:{name:'Aspartate',abbr:'Asp',letter:'D',type:'charged'},
  GAA:{name:'Glutamate',abbr:'Glu',letter:'E',type:'charged'},
  GAG:{name:'Glutamate',abbr:'Glu',letter:'E',type:'charged'},
  TGT:{name:'Cysteine',abbr:'Cys',letter:'C',type:'polar'},
  TGC:{name:'Cysteine',abbr:'Cys',letter:'C',type:'polar'},
  TGG:{name:'Tryptophan',abbr:'Trp',letter:'W',type:'nonpolar'},
  CGT:{name:'Arginine',abbr:'Arg',letter:'R',type:'charged'},
  CGC:{name:'Arginine',abbr:'Arg',letter:'R',type:'charged'},
  CGA:{name:'Arginine',abbr:'Arg',letter:'R',type:'charged'},
  CGG:{name:'Arginine',abbr:'Arg',letter:'R',type:'charged'},
  AGA:{name:'Arginine',abbr:'Arg',letter:'R',type:'charged'},
  AGG:{name:'Arginine',abbr:'Arg',letter:'R',type:'charged'},
  GGT:{name:'Glycine',abbr:'Gly',letter:'G',type:'special'},
  GGC:{name:'Glycine',abbr:'Gly',letter:'G',type:'special'},
  GGA:{name:'Glycine',abbr:'Gly',letter:'G',type:'special'},
  GGG:{name:'Glycine',abbr:'Gly',letter:'G',type:'special'},
};

/* ═══════════════════════════════════════════════════════════════════════
   BIOINFORMATICS FUNCTIONS  (mirrors the C backend logic)
═══════════════════════════════════════════════════════════════════════ */

/** Validate — only A, T, G, C allowed */
function validateDNA(seq) {
  return seq.length > 0 && /^[ATGC]+$/.test(seq);
}

/** Sanitise: strip whitespace/newlines, uppercase */
function sanitise(seq) {
  return seq.replace(/\s/g, '').toUpperCase();
}

/** DNA → mRNA  (A→U, T→A, G→C, C→G) */
function dnaToRNA(seq) {
  return seq.split('').map(b => ({A:'U',T:'A',G:'C',C:'G'}[b] || b)).join('');
}

/** Complementary DNA (A↔T, G↔C) */
function complementDNA(seq) {
  return seq.split('').map(b => ({A:'T',T:'A',G:'C',C:'G'}[b] || b)).join('');
}

/** Reverse a string */
function reverseSeq(seq) { return seq.split('').reverse().join(''); }

/** GC content percentage */
function gcContent(seq) {
  if (!seq.length) return 0;
  const gc = seq.split('').filter(b => b === 'G' || b === 'C').length;
  return (gc / seq.length * 100).toFixed(2);
}

/** Nucleotide frequency {A,T,G,C} */
function nucleotideFreq(seq) {
  const f = {A:0,T:0,G:0,C:0};
  for (const b of seq) if (b in f) f[b]++;
  return f;
}

/** Detect mutations between two equal-length sequences */
function detectMutations(ref, mut) {
  const results = [];
  const len = Math.min(ref.length, mut.length);
  for (let i = 0; i < len; i++) {
    if (ref[i] !== mut[i]) {
      results.push({ pos: i + 1, ref: ref[i], mut: mut[i] });
    }
  }
  return results;
}

/** Translate DNA coding sequence (starts at first ATG, stops at stop codon) */
function translateSequence(seq) {
  // Find start codon
  let start = -1;
  for (let i = 0; i + 2 < seq.length; i += 3) {
    if (seq.slice(i, i + 3) === 'ATG') { start = i; break; }
  }
  if (start === -1) return { error: 'No start codon (ATG) found.' };
  if ((seq.length - start) % 3 !== 0) {
    return { error: 'Sequence from ATG is not a multiple of 3.' };
  }

  const chain = [];
  for (let i = start; i + 2 < seq.length; i += 3) {
    const codon = seq.slice(i, i + 3);
    const entry = CODON_TABLE[codon] || { name:'Unknown', abbr:'???', letter:'?', type:'nonpolar' };
    chain.push({ codon, ...entry });
    if (entry.type === 'stop') break;
  }
  return { chain };
}

/** Find all occurrences of pattern in seq, returns array of start positions (1-based) */
function patternSearch(seq, pattern) {
  const positions = [];
  let idx = 0;
  while ((idx = seq.indexOf(pattern, idx)) !== -1) {
    positions.push(idx); // 0-based, converted to 1-based in render
    idx++;
  }
  return positions;
}

/* ═══════════════════════════════════════════════════════════════════════
   RENDER HELPERS
═══════════════════════════════════════════════════════════════════════ */

/** Wrap each char of seq in a <span> so mutations can be highlighted */
function renderSeqWithMutations(seq, mutPositions, cls) {
  return seq.split('').map((ch, i) => {
    if (mutPositions.has(i)) return `<span class="${cls}">${ch}</span>`;
    return ch;
  }).join('');
}

/** Highlight all pattern matches inside seq */
function renderSeqWithPattern(seq, pattern) {
  if (!pattern) return seq;
  const pLen = pattern.length;
  let html = '';
  let i = 0;
  while (i < seq.length) {
    if (seq.startsWith(pattern, i)) {
      html += `<span class="match">${seq.slice(i, i + pLen)}</span>`;
      i += pLen;
    } else {
      html += seq[i++];
    }
  }
  return html;
}

/** Build an amino acid chip element */
function aminoChip(entry) {
  const cls = entry.type === 'start' ? 'start'
            : entry.type === 'stop'  ? 'stop'
            : entry.type;
  return `<div class="amino-chip ${cls}">
    <span class="codon">${entry.codon}</span>
    <span class="abbr">${entry.abbr}</span>
    <span class="name">${entry.name}</span>
  </div>`;
}

/** Show a result-line row */
function resultRow(key, val, valClass = '') {
  return `<div class="result-line">
    <span class="result-key">${key}</span>
    <span class="result-val ${valClass} seq-display">${val}</span>
  </div>`;
}

/** Show loading skeleton */
function showLoading(el) {
  el.innerHTML = `<span class="spinner"></span>&nbsp; Analysing…`;
}

/** Animate result appearance with a fake 500 ms delay (feels responsive) */
function withDelay(fn) {
  return new Promise(resolve => setTimeout(() => { fn(); resolve(); }, 480));
}

/* ═══════════════════════════════════════════════════════════════════════
   CHART.JS — Nucleotide Frequency Donut
═══════════════════════════════════════════════════════════════════════ */
let freqChart = null;

function renderFreqChart(freq) {
  const container = document.getElementById('chart-container');
  container.classList.remove('hidden');
  const canvas = document.getElementById('freq-chart');

  if (freqChart) { freqChart.destroy(); freqChart = null; }

  freqChart = new Chart(canvas, {
    type: 'doughnut',
    data: {
      labels: ['Adenine (A)', 'Thymine (T)', 'Guanine (G)', 'Cytosine (C)'],
      datasets: [{
        data: [freq.A, freq.T, freq.G, freq.C],
        backgroundColor: ['rgba(255,107,107,.8)','rgba(255,217,61,.8)','rgba(57,255,20,.8)','rgba(77,150,255,.8)'],
        borderColor:      ['#ff6b6b','#ffd93d','#39ff14','#4d96ff'],
        borderWidth: 2,
        hoverOffset: 8,
      }]
    },
    options: {
      cutout: '65%',
      plugins: {
        legend: { position: 'bottom', labels: { color: '#7a91b5', padding: 16, font: { size: 12 } } },
        tooltip: {
          callbacks: {
            label: ctx => {
              const total = ctx.dataset.data.reduce((a,b) => a+b, 0);
              const pct = total ? (ctx.parsed / total * 100).toFixed(1) : 0;
              return ` ${ctx.parsed} bases (${pct}%)`;
            }
          }
        }
      },
      animation: { animateRotate: true, duration: 900 }
    }
  });
}

/* ═══════════════════════════════════════════════════════════════════════
   FEATURE RUNNERS
═══════════════════════════════════════════════════════════════════════ */

function getSeq() {
  return sanitise(document.getElementById('dna-input').value);
}

/* --- RNA ---------------------------------------------------------------- */
async function runRNA() {
  const box = document.getElementById('result-rna');
  const seq = getSeq();
  if (!validateDNA(seq)) { showValidationError(box); return; }
  showLoading(box);
  await withDelay(() => {
    const rna = dnaToRNA(seq);
    box.innerHTML =
      resultRow('DNA (5\'→3\')', seq) +
      resultRow('mRNA (5\'→3\')', rna, 'purple') +
      `<div class="result-line" style="margin-top:.75rem">
        <span class="result-key" style="color:#7a91b5;font-size:.75rem">Rule: A→U · T→A · G→C · C→G</span>
      </div>`;
  });
}

/* --- Complement --------------------------------------------------------- */
async function runComplement() {
  const box = document.getElementById('result-complement');
  const seq = getSeq();
  if (!validateDNA(seq)) { showValidationError(box); return; }
  showLoading(box);
  await withDelay(() => {
    const comp   = complementDNA(seq);
    const revSeq = reverseSeq(seq);
    const revComp= reverseSeq(comp);
    box.innerHTML =
      resultRow("Original   (5'→3')", seq) +
      resultRow("Complement (3'→5')", comp, 'purple') +
      resultRow("Rev. Complement (5'→3')", revComp, 'green') +
      `<div class="result-line" style="margin-top:.75rem">
        <span class="result-key" style="color:#7a91b5;font-size:.75rem">Rule: A↔T · G↔C</span>
       </div>`;
  });
}

/* --- GC Content --------------------------------------------------------- */
async function runGC() {
  const box = document.getElementById('result-gc');
  const seq = getSeq();
  if (!validateDNA(seq)) { showValidationError(box); return; }
  showLoading(box);
  await withDelay(() => {
    const gc   = parseFloat(gcContent(seq));
    const at   = (100 - gc).toFixed(2);
    const freq = nucleotideFreq(seq);
    const tot  = seq.length;

    box.innerHTML = `
      <div class="gc-stats">
        <div class="gc-stat"><span class="val gradient-text">${gc}%</span><span class="lbl">GC Content</span></div>
        <div class="gc-stat"><span class="val" style="color:var(--yellow)">${at}%</span><span class="lbl">AT Content</span></div>
        <div class="gc-stat"><span class="val" style="color:var(--txt)">${tot}</span><span class="lbl">Total bases</span></div>
      </div>
      <div class="bar-row bar-gc">
        <div class="bar-label"><span>GC</span><span>${gc}%</span></div>
        <div class="bar-track"><div class="bar-fill" style="width:0%" data-pct="${gc}"></div></div>
      </div>
      <div class="bar-row bar-at">
        <div class="bar-label"><span>AT</span><span>${at}%</span></div>
        <div class="bar-track"><div class="bar-fill" style="width:0%" data-pct="${at}"></div></div>
      </div>
      <div class="result-line" style="margin-top:.5rem">
        <span class="result-key">A</span><span class="result-val">${freq.A} &nbsp;</span>
        <span class="result-key">T</span><span class="result-val">${freq.T} &nbsp;</span>
        <span class="result-key">G</span><span class="result-val" style="color:var(--green)">${freq.G} &nbsp;</span>
        <span class="result-key">C</span><span class="result-val" style="color:var(--blue)">${freq.C}</span>
      </div>`;
    // Animate bars after render
    requestAnimationFrame(() => {
      box.querySelectorAll('.bar-fill').forEach(el => {
        el.style.width = el.dataset.pct + '%';
      });
    });
  });
}

/* --- Mutations ---------------------------------------------------------- */
async function runMutation() {
  const box = document.getElementById('result-mutation');
  const ref = sanitise(document.getElementById('ref-input').value);
  const mut = sanitise(document.getElementById('mut-input').value);

  if (!ref || !mut) { box.innerHTML = '<span class="result-placeholder">Please enter both sequences.</span>'; return; }
  if (!validateDNA(ref) || !validateDNA(mut)) {
    box.innerHTML = '<span style="color:var(--red)">⚠ Both sequences must contain only A, T, G, C.</span>';
    return;
  }
  if (ref.length !== mut.length) {
    box.innerHTML = `<span style="color:var(--red)">⚠ Length mismatch — Reference: ${ref.length} bp, Mutant: ${mut.length} bp.</span>`;
    return;
  }
  showLoading(box);
  await withDelay(() => {
    const mutations = detectMutations(ref, mut);
    const mutPositionsRef = new Set(mutations.map(m => m.pos - 1));
    const mutPositionsMut = new Set(mutations.map(m => m.pos - 1));

    const refHtml = renderSeqWithMutations(ref, mutPositionsRef, 'mut-ref');
    const mutHtml = renderSeqWithMutations(mut, mutPositionsMut, 'mut-alt');

    let tableRows = mutations.length === 0
      ? '<tr><td colspan="4" style="text-align:center;color:var(--green);padding:1rem">✓ Sequences are identical — no mutations found</td></tr>'
      : mutations.map(m => `<tr>
          <td><span class="badge-pos">${m.pos}</span></td>
          <td><span class="badge-ref">${m.ref}</span></td>
          <td><span class="badge-alt">${m.mut}</span></td>
          <td><span class="badge-type">SNP</span></td>
        </tr>`).join('');

    box.innerHTML = `
      <div class="result-line"><span class="result-key">Reference</span><span class="result-val seq-display">${refHtml}</span></div>
      <div class="result-line"><span class="result-key">Mutant</span><span class="result-val seq-display" style="color:var(--orange)">${mutHtml}</span></div>
      <table class="mut-table">
        <thead><tr><th>Position</th><th>Ref Base</th><th>Mut Base</th><th>Type</th></tr></thead>
        <tbody>${tableRows}</tbody>
      </table>
      <div class="result-line" style="margin-top:.75rem">
        <span class="result-key">Total mutations</span>
        <span class="result-val" style="color:${mutations.length > 0 ? 'var(--red)' : 'var(--green)'}">${mutations.length}</span>
      </div>`;
  });
}

/* --- Translation ------------------------------------------------------- */
async function runTranslate() {
  const box = document.getElementById('result-translate');
  const seq = getSeq();
  if (!validateDNA(seq)) { showValidationError(box); return; }
  showLoading(box);
  await withDelay(() => {
    const result = translateSequence(seq);
    if (result.error) {
      box.innerHTML = `<span style="color:var(--red)">⚠ ${result.error}</span>`;
      return;
    }
    const chips = result.chain.map((entry, i) => {
      const arrow = (i < result.chain.length - 1) ? '<span class="amino-arrow">→</span>' : '';
      return aminoChip(entry) + arrow;
    }).join('');

    const oneLetter = result.chain.map(e => e.letter).join('');
    box.innerHTML = `
      <div class="amino-chain">${chips}</div>
      <div class="result-line" style="margin-top:1rem">
        <span class="result-key">Single-letter</span>
        <span class="result-val" style="letter-spacing:.08em">${oneLetter}</span>
      </div>
      <div class="result-line">
        <span class="result-key">Chain length</span>
        <span class="result-val">${result.chain.length} amino acids</span>
      </div>`;
  });
}

/* --- Pattern Search ----------------------------------------------------- */
async function runPattern() {
  const box     = document.getElementById('result-pattern');
  const seq     = getSeq();
  const pattern = sanitise(document.getElementById('pattern-input').value);

  if (!validateDNA(seq)) { showValidationError(box); return; }
  if (!pattern) { box.innerHTML = '<span class="result-placeholder">Enter a pattern above.</span>'; return; }
  if (!validateDNA(pattern)) {
    box.innerHTML = '<span style="color:var(--red)">⚠ Pattern must contain only A, T, G, C.</span>';
    return;
  }
  showLoading(box);
  await withDelay(() => {
    const positions = patternSearch(seq, pattern);
    const highlighted = renderSeqWithPattern(seq, pattern);

    const posText = positions.length > 0
      ? positions.map(p => `<span class="badge-pos">${p + 1}</span>`).join(' ')
      : '<span style="color:var(--red)">Not found</span>';

    box.innerHTML = `
      <div class="result-line"><span class="result-key">Sequence</span></div>
      <div class="seq-display" style="margin-bottom:1rem;color:var(--txt)">${highlighted}</div>
      <div class="result-line"><span class="result-key">Match count</span>
        <span class="result-val" style="color:${positions.length > 0 ? 'var(--green)' : 'var(--red)'}">${positions.length}</span>
      </div>
      <div class="result-line"><span class="result-key">Positions (1-based)</span></div>
      <div style="display:flex;flex-wrap:wrap;gap:.4rem;margin-top:.25rem">${posText}</div>`;
  });
}

/* --- Reverse ------------------------------------------------------------ */
async function runReverse() {
  const box = document.getElementById('result-reverse');
  const seq = getSeq();
  if (!validateDNA(seq)) { showValidationError(box); return; }
  showLoading(box);
  await withDelay(() => {
    const rev     = reverseSeq(seq);
    const comp    = complementDNA(seq);
    const revComp = reverseSeq(comp);
    box.innerHTML =
      resultRow("Original             ", seq) +
      resultRow("Reversed (3'→5')     ", rev, 'purple') +
      resultRow("Rev. Complement (5'→3')", revComp, 'green');
  });
}

/* --- Frequency ---------------------------------------------------------- */
async function runFrequency() {
  const box = document.getElementById('result-frequency');
  const seq = getSeq();
  if (!validateDNA(seq)) { showValidationError(box); return; }
  showLoading(box);
  await withDelay(() => {
    const freq  = nucleotideFreq(seq);
    const total = seq.length;
    const pct   = base => total ? (freq[base] / total * 100).toFixed(1) : 0;
    const maxF  = Math.max(freq.A, freq.T, freq.G, freq.C) || 1;

    const barsHtml = ['A','T','G','C'].map(base => `
      <div class="freq-col" data-base="${base}">
        <div class="freq-pct">${pct(base)}%</div>
        <div class="freq-bar-wrap">
          <div class="freq-bar" style="height:0px" data-h="${(freq[base]/maxF*120).toFixed(0)}"></div>
        </div>
        <div class="freq-base">${base}</div>
        <div class="freq-count">${freq[base]}</div>
      </div>`).join('');

    box.innerHTML = `<div class="freq-bars">${barsHtml}</div>`;

    // Animate bars
    requestAnimationFrame(() => {
      box.querySelectorAll('.freq-bar').forEach(el => {
        el.style.height = el.dataset.h + 'px';
      });
    });

    renderFreqChart(freq);
  });
}

/* --- File --------------------------------------------------------------- */
function handleFileContent(text, filename) {
  const box = document.getElementById('result-file');
  box.classList.remove('hidden');

  // Parse FASTA-like: skip # and > lines
  const lines = text.split('\n')
    .map(l => l.trim())
    .filter(l => l && !l.startsWith('#') && !l.startsWith('>'));

  if (lines.length === 0) {
    box.innerHTML = '<span style="color:var(--red)">⚠ No valid sequences found in file.</span>';
    return;
  }

  let html = `<div style="margin-bottom:.75rem;color:var(--txt2);font-size:.8rem">📄 ${filename} — ${lines.length} sequence(s) found</div>`;

  lines.forEach((line, idx) => {
    const seq = sanitise(line);
    const valid = validateDNA(seq);
    const gc    = valid ? gcContent(seq) : '—';
    const freq  = valid ? nucleotideFreq(seq) : null;

    html += `<div style="margin-bottom:1rem;padding:.75rem;background:rgba(0,0,0,.2);border-radius:8px;border:1px solid var(--border)">
      <div style="font-size:.72rem;color:var(--txt2);margin-bottom:.4rem;text-transform:uppercase;letter-spacing:.06em">Sequence #${idx + 1}</div>
      <div class="seq-display" style="color:${valid ? 'var(--cyan)' : 'var(--red)'};margin-bottom:.5rem">${seq}</div>
      ${valid ? `<div style="font-size:.78rem;color:var(--txt2)">Length: ${seq.length} bp &nbsp;|&nbsp; GC: ${gc}% &nbsp;|&nbsp; A:${freq.A} T:${freq.T} G:${freq.G} C:${freq.C}</div>` : '<div style="color:var(--red);font-size:.78rem">⚠ Invalid DNA — contains non-ATGC characters</div>'}
    </div>`;

    // Auto-load first valid sequence into main input
    if (idx === 0 && valid) {
      document.getElementById('dna-input').value = seq;
      updateLiveStats();
    }
  });

  box.innerHTML = html;
}

/* ── Validation error helper ──────────────────────────────────────── */
function showValidationError(box) {
  box.innerHTML = '<span style="color:var(--red)">⚠ Invalid or empty DNA sequence. Only A, T, G, C are allowed.</span>';
}

/* ═══════════════════════════════════════════════════════════════════════
   CANVAS ANIMATIONS
═══════════════════════════════════════════════════════════════════════ */

/* ── DNA Helix (hero) ─────────────────────────────────────────────── */
function initHelixCanvas() {
  const canvas = document.getElementById('helix-canvas');
  if (!canvas) return;
  const ctx = canvas.getContext('2d');
  let t = 0;
  const pairs = [['A','T'],['T','A'],['G','C'],['C','G'],['A','T'],['G','C'],['T','A'],['C','G']];
  const colors = {A:'#ff6b6b',T:'#ffd93d',G:'#39ff14',C:'#4d96ff'};

  function resize() { canvas.width = canvas.offsetWidth; canvas.height = canvas.offsetHeight; }
  window.addEventListener('resize', resize); resize();

  function draw() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    const cx  = canvas.width / 2;
    const amp = Math.min(canvas.width * 0.32, 75);
    const sp  = 32;
    const num = Math.ceil(canvas.height / sp) + 2;

    // Backbone paths
    ctx.beginPath();
    for (let i = 0; i <= num; i++) {
      const y = i * sp - (t % sp);
      const x = cx + Math.sin((i * sp + t) * 0.045) * amp;
      i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
    }
    ctx.strokeStyle = 'rgba(0,245,212,0.45)'; ctx.lineWidth = 2;
    ctx.shadowColor = 'rgba(0,245,212,0.6)'; ctx.shadowBlur = 6; ctx.stroke(); ctx.shadowBlur = 0;

    ctx.beginPath();
    for (let i = 0; i <= num; i++) {
      const y = i * sp - (t % sp);
      const x = cx - Math.sin((i * sp + t) * 0.045) * amp;
      i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
    }
    ctx.strokeStyle = 'rgba(139,92,246,0.45)'; ctx.lineWidth = 2;
    ctx.shadowColor = 'rgba(139,92,246,0.6)'; ctx.shadowBlur = 6; ctx.stroke(); ctx.shadowBlur = 0;

    // Rungs & nucleotide dots
    for (let i = 0; i < num; i++) {
      const y     = i * sp - (t % sp);
      const phase = (i * sp + t) * 0.045;
      const x1    = cx + Math.sin(phase) * amp;
      const x2    = cx - Math.sin(phase) * amp;
      const op    = 0.25 + 0.5 * Math.abs(Math.sin(phase));
      const pairIdx = Math.abs(i) % pairs.length;
      const [lB, rB] = pairs[pairIdx];

      // Rung
      const grad = ctx.createLinearGradient(x1, y, x2, y);
      grad.addColorStop(0, `rgba(0,245,212,${op})`);
      grad.addColorStop(1, `rgba(139,92,246,${op})`);
      ctx.beginPath(); ctx.moveTo(x1, y); ctx.lineTo(x2, y);
      ctx.strokeStyle = grad; ctx.lineWidth = 1.5;
      ctx.globalAlpha = op + 0.2; ctx.stroke(); ctx.globalAlpha = 1;

      // Left dot
      ctx.shadowColor = colors[lB]; ctx.shadowBlur = 8;
      ctx.beginPath(); ctx.arc(x1, y, 5.5, 0, Math.PI * 2);
      ctx.fillStyle = colors[lB]; ctx.globalAlpha = op + 0.3; ctx.fill();

      // Right dot
      ctx.beginPath(); ctx.arc(x2, y, 5.5, 0, Math.PI * 2);
      ctx.fillStyle = colors[rB]; ctx.fill();
      ctx.shadowBlur = 0; ctx.globalAlpha = 1;

      // Base labels
      ctx.font = 'bold 7px Inter,sans-serif'; ctx.fillStyle = '#000';
      ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
      ctx.globalAlpha = op + 0.4;
      ctx.fillText(lB, x1, y); ctx.fillText(rB, x2, y);
      ctx.globalAlpha = 1;
    }
    t += 0.55;
    requestAnimationFrame(draw);
  }
  draw();
}

/* ── Particle Field (background) ──────────────────────────────────── */
function initParticles() {
  const canvas = document.getElementById('particle-canvas');
  if (!canvas) return;
  const ctx = canvas.getContext('2d');

  function resize() { canvas.width = window.innerWidth; canvas.height = window.innerHeight; }
  window.addEventListener('resize', resize); resize();

  const bases   = ['A','T','G','C'];
  const bColors = ['rgba(255,107,107,0.35)','rgba(255,217,61,0.35)','rgba(57,255,20,0.35)','rgba(77,150,255,0.35)'];
  const count   = Math.min(60, Math.floor(window.innerWidth / 22));

  const particles = Array.from({length: count}, () => ({
    x: Math.random() * canvas.width,
    y: Math.random() * canvas.height,
    s: 10 + Math.random() * 8,
    sp: 0.2 + Math.random() * 0.5,
    op: 0.15 + Math.random() * 0.35,
    b: Math.floor(Math.random() * 4),
  }));

  function draw() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    particles.forEach(p => {
      ctx.font = `${p.s}px 'JetBrains Mono',monospace`;
      ctx.fillStyle = bColors[p.b];
      ctx.globalAlpha = p.op;
      ctx.fillText(bases[p.b], p.x, p.y);
      p.y -= p.sp;
      if (p.y < -20) { p.y = canvas.height + 20; p.x = Math.random() * canvas.width; }
    });
    ctx.globalAlpha = 1;
    requestAnimationFrame(draw);
  }
  draw();
}

/* ═══════════════════════════════════════════════════════════════════════
   LIVE STATS  (updates as user types)
═══════════════════════════════════════════════════════════════════════ */
function updateLiveStats() {
  const raw  = document.getElementById('dna-input').value;
  const seq  = sanitise(raw);
  const lenEl = document.getElementById('live-length');
  const gcEl  = document.getElementById('live-gc');

  lenEl.textContent = `${seq.length} bp`;

  if (seq.length === 0) {
    gcEl.textContent = 'GC: —';
  } else if (validateDNA(seq)) {
    gcEl.textContent = `GC: ${gcContent(seq)}%`;
    gcEl.style.color = '';
  } else {
    gcEl.textContent = 'Invalid';
    gcEl.style.color = 'var(--red)';
  }
}

/* ═══════════════════════════════════════════════════════════════════════
   TAB SWITCHING
═══════════════════════════════════════════════════════════════════════ */
function initTabs() {
  document.querySelectorAll('.tab').forEach(btn => {
    btn.addEventListener('click', () => {
      document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
      document.querySelectorAll('.panel').forEach(p => p.classList.remove('active'));
      btn.classList.add('active');
      document.getElementById(`panel-${btn.dataset.tab}`).classList.add('active');
    });
  });
}

/* ═══════════════════════════════════════════════════════════════════════
   VALIDATION BUTTON
═══════════════════════════════════════════════════════════════════════ */
function initValidation() {
  document.getElementById('validate-btn').addEventListener('click', () => {
    const seq = getSeq();
    const msg = document.getElementById('validation-msg');
    msg.classList.remove('hidden', 'valid', 'invalid');
    if (!seq) {
      msg.textContent = '⚠ Sequence is empty. Please enter a DNA sequence.';
      msg.classList.add('invalid');
    } else if (validateDNA(seq)) {
      msg.textContent = `✓ Valid DNA sequence — ${seq.length} bp, GC: ${gcContent(seq)}%`;
      msg.classList.add('valid');
    } else {
      const bad = [...new Set(seq.split('').filter(c => !'ATGC'.includes(c)))].join(', ');
      msg.textContent = `✗ Invalid characters detected: ${bad}. Only A, T, G, C are allowed.`;
      msg.classList.add('invalid');
    }
  });

  document.getElementById('clear-btn').addEventListener('click', () => {
    document.getElementById('dna-input').value = '';
    document.getElementById('validation-msg').classList.add('hidden');
    updateLiveStats();
  });

  document.getElementById('sample-btn').addEventListener('click', () => {
    document.getElementById('dna-input').value = 'ATGTTTTGGAAACGCCCGGCTGAATTTAAGCATAAGGTGAGTTTCGAGTAA';
    document.getElementById('validation-msg').classList.add('hidden');
    updateLiveStats();
  });

  document.getElementById('dna-input').addEventListener('input', updateLiveStats);
}

/* ═══════════════════════════════════════════════════════════════════════
   RUN BUTTONS  (generic + per-panel)
═══════════════════════════════════════════════════════════════════════ */
function initRunButtons() {
  const actions = {
    runRNA:       runRNA,
    runComplement:runComplement,
    runGC:        runGC,
    runMutation:  runMutation,
    runTranslate: runTranslate,
    runReverse:   runReverse,
    runFrequency: runFrequency,
  };
  document.querySelectorAll('.run-btn[data-action]').forEach(btn => {
    btn.addEventListener('click', () => actions[btn.dataset.action]?.());
  });
  document.getElementById('run-pattern-btn').addEventListener('click', runPattern);
  document.getElementById('pattern-input').addEventListener('keydown', e => {
    if (e.key === 'Enter') runPattern();
  });
}

/* ═══════════════════════════════════════════════════════════════════════
   FILE HANDLING
═══════════════════════════════════════════════════════════════════════ */
function initFileHandlers() {
  // Top-level file button (loads to textarea)
  document.getElementById('file-upload').addEventListener('change', e => {
    const file = e.target.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = ev => {
      const lines = ev.target.result.split('\n')
        .map(l => l.trim()).filter(l => l && !l.startsWith('#') && !l.startsWith('>'));
      if (lines[0]) {
        document.getElementById('dna-input').value = sanitise(lines[0]);
        updateLiveStats();
      }
    };
    reader.readAsText(file);
    e.target.value = '';
  });

  // File panel uploader
  document.getElementById('file-upload-panel').addEventListener('change', e => {
    const file = e.target.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = ev => handleFileContent(ev.target.result, file.name);
    reader.readAsText(file);
    e.target.value = '';
  });

  // Drag & drop zone
  const dropZone = document.getElementById('drop-zone');
  dropZone.addEventListener('dragover', e => { e.preventDefault(); dropZone.classList.add('drag-over'); });
  dropZone.addEventListener('dragleave', () => dropZone.classList.remove('drag-over'));
  dropZone.addEventListener('drop', e => {
    e.preventDefault(); dropZone.classList.remove('drag-over');
    const file = e.dataTransfer.files[0];
    if (!file || !file.name.endsWith('.txt')) {
      alert('Please drop a .txt file.');
      return;
    }
    const reader = new FileReader();
    reader.onload = ev => handleFileContent(ev.target.result, file.name);
    reader.readAsText(file);
  });
}

/* ═══════════════════════════════════════════════════════════════════════
   NAVBAR SCROLL EFFECT & HAMBURGER
═══════════════════════════════════════════════════════════════════════ */
function initNavbar() {
  const nav = document.getElementById('navbar');
  window.addEventListener('scroll', () => {
    nav.classList.toggle('scrolled', window.scrollY > 40);
  });
  document.getElementById('hamburger').addEventListener('click', () => {
    document.querySelector('.nav-links').classList.toggle('open');
  });
}

/* ═══════════════════════════════════════════════════════════════════════
   INITIALISE
═══════════════════════════════════════════════════════════════════════ */
document.addEventListener('DOMContentLoaded', () => {
  initParticles();
  initHelixCanvas();
  initNavbar();
  initTabs();
  initValidation();
  initRunButtons();
  initFileHandlers();
  updateLiveStats();
});
