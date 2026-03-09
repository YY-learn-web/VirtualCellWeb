import React, { useMemo, useState } from 'react';
import { Dna, Network, Beaker, KeyRound, Download, FlaskConical } from 'lucide-react';
import { ResponsiveContainer, BarChart, Bar, XAxis, YAxis, Tooltip, CartesianGrid } from 'recharts';
import { PRESET_CELL_LINES, MAX_TIME_HOURS, MAX_DOSE_UM } from '../constants';
import { Button } from './Button';
import { Input, RangeSlider } from './Input';
import { ResultsChart } from './ResultsChart';
import { ASCENDService } from '../services/apiLive';
import { buildApiUrl } from '../services/apiConfig';
import { EnrichmentTerm, GeneExpressionResult, IC50Prediction, KeyGeneScore, SingleMetabolomicsResponse, StringNetworkResponse } from '../types';

interface SingleWorkflowPageProps {
  isDark: boolean;
}

export const SingleWorkflowPage: React.FC<SingleWorkflowPageProps> = ({ isDark }) => {
  const [smiles, setSmiles] = useState<string>('CC(=O)OC1=CC=CC=C1C(=O)O');
  const [time, setTime] = useState<number>(24);
  const [dose, setDose] = useState<number>(10);
  const [cellLineId, setCellLineId] = useState<string>(PRESET_CELL_LINES[0].id);
  const [enrichmentLibrary, setEnrichmentLibrary] = useState<string>('GO_Biological_Process_2021');
  const [topGenesForEnrichment, setTopGenesForEnrichment] = useState<number>(200);

  const [transcriptomics, setTranscriptomics] = useState<GeneExpressionResult[] | null>(null);
  const [enrichment, setEnrichment] = useState<EnrichmentTerm[]>([]);
  const [ppi, setPpi] = useState<StringNetworkResponse | null>(null);
  const [metabolomics, setMetabolomics] = useState<SingleMetabolomicsResponse | null>(null);
  const [ic50, setIc50] = useState<IC50Prediction | null>(null);
  const [keyGenes, setKeyGenes] = useState<KeyGeneScore[]>([]);

  const [loadingStep, setLoadingStep] = useState<'step1' | 'step2' | 'step3' | 'step4' | 'step4_ic50' | 'step5' | null>(null);
  const [errorMsg, setErrorMsg] = useState<string | null>(null);

  const rankedGenes = useMemo(() => {
    if (!transcriptomics?.length) return [];
    return [...transcriptomics].sort((a, b) => Math.abs(b.expressionLevel) - Math.abs(a.expressionLevel));
  }, [transcriptomics]);

  const enrichmentChart = useMemo(() => {
    return enrichment
      .slice(0, 12)
      .map((term) => ({
        term: term.term,
        score: -Math.log10((term.adjP ?? term.pValue) + 1e-12),
      }))
      .reverse();
  }, [enrichment]);

  const topFluxEntries = useMemo(() => {
    if (!metabolomics?.flux) return [];
    return Object.entries(metabolomics.flux).slice(0, 20);
  }, [metabolomics]);

  const downloadUrl = (path: string) => buildApiUrl(path);

  const runStep1Transcriptomics = async () => {
    setLoadingStep('step1');
    setErrorMsg(null);
    setEnrichment([]);
    setPpi(null);
    setMetabolomics(null);
    setIc50(null);
    setKeyGenes([]);
    try {
      const response = await ASCENDService.predictExpression({
        smiles,
        time,
        dose,
        cellLineMode: 'preset',
        cellLineId,
      });
      setTranscriptomics(response.data);
    } catch (err: any) {
      setTranscriptomics(null);
      setErrorMsg(err?.message ?? 'Step 1 failed.');
    } finally {
      setLoadingStep(null);
    }
  };

  const runStep2Enrichment = async () => {
    if (!rankedGenes.length) {
      setErrorMsg('Run Step 1 first.');
      return;
    }
    setLoadingStep('step2');
    setErrorMsg(null);
    try {
      const geneIds = rankedGenes.slice(0, topGenesForEnrichment).map((g) => g.geneId);
      const response = await ASCENDService.runEnrichment({
        geneIds,
        library: enrichmentLibrary,
        topN: 20,
      });
      setEnrichment(response.data);
    } catch (err: any) {
      setEnrichment([]);
      setErrorMsg(err?.message ?? 'Step 2 failed.');
    } finally {
      setLoadingStep(null);
    }
  };

  const runStep3Ppi = async () => {
    if (!transcriptomics?.length) {
      setErrorMsg('Run Step 1 first.');
      return;
    }
    setLoadingStep('step3');
    setErrorMsg(null);
    try {
      const response = await ASCENDService.fetchStringNetwork({
        genes: transcriptomics,
        minExpression: 1.5,
        species: 9606,
        requiredScore: 400,
        networkType: 'functional',
        maxGenes: 200,
      });
      setPpi(response);
    } catch (err: any) {
      setPpi(null);
      setErrorMsg(err?.message ?? 'Step 3 failed.');
    } finally {
      setLoadingStep(null);
    }
  };

  const runStep4Metabolomics = async () => {
    if (!transcriptomics?.length) {
      setErrorMsg('Run Step 1 first.');
      return;
    }
    setLoadingStep('step4');
    setErrorMsg(null);
    try {
      const response = await ASCENDService.runSingleWorkflowMetabolomics({
        expression: transcriptomics,
        epochs: 100,
        imputation: false,
      });
      setMetabolomics(response);
    } catch (err: any) {
      setMetabolomics(null);
      setErrorMsg(err?.message ?? 'Step 4 failed.');
    } finally {
      setLoadingStep(null);
    }
  };

  const runStep4Ic50 = async () => {
    if (!transcriptomics?.length) {
      setErrorMsg('Run Step 1 first.');
      return;
    }
    setLoadingStep('step4_ic50');
    setErrorMsg(null);
    try {
      const response = await ASCENDService.predictIC50({
        cellLineId,
        smiles,
        expression: transcriptomics.map(({ geneId, expressionLevel }) => ({ geneId, expressionLevel })),
        metadata: {
          timeHours: time,
          doseUm: dose,
        },
      });
      setIc50(response.data);
    } catch (err: any) {
      setIc50(null);
      setErrorMsg(err?.message ?? 'IC50 prediction failed.');
    } finally {
      setLoadingStep(null);
    }
  };

  const runStep5KeyGenes = async () => {
    if (!transcriptomics?.length) {
      setErrorMsg('Run Step 1 first.');
      return;
    }
    setLoadingStep('step5');
    setErrorMsg(null);
    try {
      const response = await ASCENDService.identifyKeyGenes({
        expression: transcriptomics,
        enrichment,
        ppiGenes: ppi?.genesUsed ?? [],
        metabolomicsFlux: metabolomics?.flux ?? {},
        topN: 5,
      });
      setKeyGenes(response.topGenes);
    } catch (err: any) {
      setKeyGenes([]);
      setErrorMsg(err?.message ?? 'Step 5 failed.');
    } finally {
      setLoadingStep(null);
    }
  };

  const axisColor = isDark ? '#94a3b8' : '#64748b';
  const gridColor = isDark ? '#1e293b' : '#e2e8f0';

  return (
    <div className="space-y-8">
      <div className="text-center max-w-4xl mx-auto">
        <h1 className="text-4xl md:text-5xl font-bold text-slate-900 dark:text-white mb-3 tracking-tight">
          Single-Sample <span className="text-transparent bg-clip-text bg-gradient-to-r from-cyan-600 to-blue-600 dark:from-cyan-400 dark:to-blue-500">Integrated Workflow</span>
        </h1>
        <p className="text-sm md:text-base text-slate-500 dark:text-slate-400">
          One ordered workflow: transcriptomics, enrichment, PPI network, metabolomics quantification, and key-gene ranking.
        </p>
      </div>

      {errorMsg && (
        <div className="max-w-5xl mx-auto p-3 rounded-lg bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-500/40 text-sm text-red-600 dark:text-red-300">
          {errorMsg}
        </div>
      )}

      <div className="max-w-5xl mx-auto space-y-6">
        <section className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md">
          <div className="flex items-center justify-between gap-4 mb-4">
            <h2 className="text-lg font-bold text-slate-900 dark:text-white flex items-center gap-2">
              <Dna className="w-5 h-5 text-cyan-600 dark:text-cyan-400" /> Step 1: Transcriptomics Prediction
            </h2>
            <Button onClick={runStep1Transcriptomics} isLoading={loadingStep === 'step1'}>
              Run Step 1
            </Button>
          </div>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <Input label="Drug SMILES" value={smiles} onChange={(e) => setSmiles(e.target.value)} />
            <div>
              <label className="block mb-2 text-xs font-mono tracking-widest text-slate-500 dark:text-cyan-500 uppercase font-bold">
                Cell Line
              </label>
              <select
                value={cellLineId}
                onChange={(e) => setCellLineId(e.target.value)}
                className="w-full px-4 py-3 text-sm text-slate-900 dark:text-white bg-slate-50 dark:bg-slate-950/50 border border-slate-200 dark:border-slate-700 rounded-lg"
              >
                {PRESET_CELL_LINES.map((c) => (
                  <option key={c.id} value={c.id}>{c.name}</option>
                ))}
              </select>
            </div>
          </div>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6 mt-4">
            <RangeSlider label="Time" min={0} max={MAX_TIME_HOURS} unit="h" value={time} onChange={(e) => setTime(Number(e.target.value))} />
            <RangeSlider label="Dose" min={0} max={MAX_DOSE_UM} unit="uM" value={dose} onChange={(e) => setDose(Number(e.target.value))} />
          </div>
          {transcriptomics?.length ? (
            <div className="mt-4">
              <ResultsChart data={transcriptomics} isDark={isDark} />
            </div>
          ) : null}
        </section>

        <section className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md">
          <div className="flex items-center justify-between gap-4 mb-4">
            <h2 className="text-lg font-bold text-slate-900 dark:text-white flex items-center gap-2">
              <FlaskConical className="w-5 h-5 text-cyan-600 dark:text-cyan-400" /> Step 2: Gene Enrichment
            </h2>
            <Button onClick={runStep2Enrichment} isLoading={loadingStep === 'step2'}>
              Run Step 2
            </Button>
          </div>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-4">
            <Input
              label="Top Differential Genes"
              type="number"
              min={20}
              max={1000}
              value={topGenesForEnrichment}
              onChange={(e) => setTopGenesForEnrichment(Number(e.target.value))}
            />
            <div>
              <label className="block mb-2 text-xs font-mono tracking-widest text-slate-500 dark:text-cyan-500 uppercase font-bold">
                GO Library
              </label>
              <select
                value={enrichmentLibrary}
                onChange={(e) => setEnrichmentLibrary(e.target.value)}
                className="w-full px-4 py-3 text-sm text-slate-900 dark:text-white bg-slate-50 dark:bg-slate-950/50 border border-slate-200 dark:border-slate-700 rounded-lg"
              >
                <option value="GO_Biological_Process_2021">GO Biological Process</option>
                <option value="GO_Molecular_Function_2021">GO Molecular Function</option>
                <option value="GO_Cellular_Component_2021">GO Cellular Component</option>
              </select>
            </div>
          </div>
          {enrichment.length > 0 ? (
            <div className="h-[320px]">
              <ResponsiveContainer width="100%" height="100%">
                <BarChart data={enrichmentChart} layout="vertical" margin={{ left: 12, right: 24 }}>
                  <CartesianGrid stroke={gridColor} strokeDasharray="3 3" />
                  <XAxis type="number" stroke={axisColor} />
                  <YAxis dataKey="term" type="category" width={180} stroke={axisColor} />
                  <Tooltip />
                  <Bar dataKey="score" fill="#06b6d4" />
                </BarChart>
              </ResponsiveContainer>
            </div>
          ) : null}
        </section>

        <section className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md">
          <div className="flex items-center justify-between gap-4 mb-4">
            <h2 className="text-lg font-bold text-slate-900 dark:text-white flex items-center gap-2">
              <Network className="w-5 h-5 text-cyan-600 dark:text-cyan-400" /> Step 3: Protein-Protein Interaction
            </h2>
            <Button onClick={runStep3Ppi} isLoading={loadingStep === 'step3'}>
              Run Step 3
            </Button>
          </div>
          {ppi ? (
            <div className="space-y-3">
              <p className="text-xs text-slate-500 dark:text-slate-400">
                Filtered genes (&gt;1.5): {ppi.filteredCount} | Mapped genes: {ppi.mappedCount}
              </p>
              {ppi.link ? (
                <a href={ppi.link} target="_blank" rel="noreferrer" className="text-xs text-cyan-600 dark:text-cyan-400 underline">
                  Open interactive STRING view
                </a>
              ) : null}
              <div className="w-full overflow-x-auto rounded-lg border border-slate-200 dark:border-slate-800 bg-slate-50 dark:bg-slate-950/40 p-2">
                <img src={ppi.image} alt="STRING network" className="w-full h-auto" />
              </div>
            </div>
          ) : null}
        </section>

        <section className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md">
          <div className="flex items-center justify-between gap-4 mb-4">
            <h2 className="text-lg font-bold text-slate-900 dark:text-white flex items-center gap-2">
              <Beaker className="w-5 h-5 text-cyan-600 dark:text-cyan-400" /> Step 4: Metabolomics Quantification
            </h2>
            <Button onClick={runStep4Metabolomics} isLoading={loadingStep === 'step4'}>
              Run Step 4
            </Button>
          </div>
          {metabolomics ? (
            <div className="space-y-4">
              <div className="flex flex-wrap gap-3">
                <a href={downloadUrl(metabolomics.files.flux)} target="_blank" rel="noreferrer">
                  <Button variant="secondary" className="flex items-center gap-2">
                    <Download className="w-4 h-4" /> Download Flux CSV
                  </Button>
                </a>
                <a href={downloadUrl(metabolomics.files.balance)} target="_blank" rel="noreferrer">
                  <Button variant="secondary" className="flex items-center gap-2">
                    <Download className="w-4 h-4" /> Download Balance CSV
                  </Button>
                </a>
              </div>
              <div className="overflow-x-auto border border-slate-200 dark:border-slate-800 rounded-lg">
                <table className="min-w-full text-left text-sm">
                  <thead className="bg-slate-50 dark:bg-slate-800 text-slate-500 dark:text-slate-400 uppercase text-xs">
                    <tr>
                      <th className="px-4 py-2">Module</th>
                      <th className="px-4 py-2">Flux</th>
                    </tr>
                  </thead>
                  <tbody className="text-slate-700 dark:text-slate-300">
                    {topFluxEntries.map(([module, value]) => (
                      <tr key={module} className="border-t border-slate-200 dark:border-slate-800">
                        <td className="px-4 py-2">{module}</td>
                        <td className="px-4 py-2">{value.toFixed(6)}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          ) : null}
        </section>

        <section className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md">
          <div className="flex items-center justify-between gap-4 mb-4">
            <h2 className="text-lg font-bold text-slate-900 dark:text-white flex items-center gap-2">
              <FlaskConical className="w-5 h-5 text-cyan-600 dark:text-cyan-400" /> Step 4.5: IC50 Prediction
            </h2>
            <Button onClick={runStep4Ic50} isLoading={loadingStep === 'step4_ic50'}>
              Run IC50
            </Button>
          </div>
          {ic50 ? (
            <div className="grid grid-cols-1 sm:grid-cols-3 gap-4">
              <div className="p-4 rounded-lg bg-gradient-to-r from-cyan-500 to-blue-500 text-white">
                <p className="text-xs uppercase tracking-wider opacity-90">IC50</p>
                <p className="text-2xl font-bold">{ic50.ic50.toFixed(4)} {ic50.unit}</p>
              </div>
              <div className="p-4 rounded-lg bg-slate-100 dark:bg-slate-800">
                <p className="text-xs uppercase tracking-wider text-slate-500 dark:text-slate-400">ln(IC50)</p>
                <p className="text-2xl font-bold text-slate-900 dark:text-white">{ic50.lnIc50.toFixed(4)}</p>
              </div>
              <div className="p-4 rounded-lg bg-slate-100 dark:bg-slate-800">
                <p className="text-xs uppercase tracking-wider text-slate-500 dark:text-slate-400">Confidence</p>
                <p className="text-2xl font-bold text-slate-900 dark:text-white">
                  {ic50.confidence !== undefined ? `${(ic50.confidence * 100).toFixed(1)}%` : 'N/A'}
                </p>
              </div>
            </div>
          ) : (
            <p className="text-sm text-slate-500 dark:text-slate-400">
              Uses transcriptomics output from Step 1 for single-sample IC50 estimation.
            </p>
          )}
        </section>

        <section className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md">
          <div className="flex items-center justify-between gap-4 mb-4">
            <h2 className="text-lg font-bold text-slate-900 dark:text-white flex items-center gap-2">
              <KeyRound className="w-5 h-5 text-cyan-600 dark:text-cyan-400" /> Step 5: Key Gene Identification
            </h2>
            <Button onClick={runStep5KeyGenes} isLoading={loadingStep === 'step5'}>
              Run Step 5
            </Button>
          </div>
          {keyGenes.length > 0 ? (
            <div className="overflow-x-auto border border-slate-200 dark:border-slate-800 rounded-lg">
              <table className="min-w-full text-left text-sm">
                <thead className="bg-slate-50 dark:bg-slate-800 text-slate-500 dark:text-slate-400 uppercase text-xs">
                  <tr>
                    <th className="px-4 py-2">Rank</th>
                    <th className="px-4 py-2">Gene</th>
                    <th className="px-4 py-2">Final Score</th>
                    <th className="px-4 py-2">Expression</th>
                    <th className="px-4 py-2">Enrichment</th>
                    <th className="px-4 py-2">PPI</th>
                    <th className="px-4 py-2">Metabolomics</th>
                  </tr>
                </thead>
                <tbody className="text-slate-700 dark:text-slate-300">
                  {keyGenes.map((row, idx) => (
                    <tr key={row.gene} className="border-t border-slate-200 dark:border-slate-800">
                      <td className="px-4 py-2">{idx + 1}</td>
                      <td className="px-4 py-2 font-semibold">{row.gene}</td>
                      <td className="px-4 py-2">{row.score.toFixed(4)}</td>
                      <td className="px-4 py-2">{row.expressionScore.toFixed(3)}</td>
                      <td className="px-4 py-2">{row.enrichmentScore.toFixed(3)}</td>
                      <td className="px-4 py-2">{row.ppiScore.toFixed(3)}</td>
                      <td className="px-4 py-2">{row.metabolomicsScore.toFixed(3)}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          ) : null}
        </section>
      </div>
    </div>
  );
};
