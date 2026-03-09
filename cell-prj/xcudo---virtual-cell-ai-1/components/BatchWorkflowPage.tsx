import React, { useMemo, useState } from 'react';
import { Activity, Download, FlaskConical, Network, Beaker, KeyRound } from 'lucide-react';
import { ResponsiveContainer, BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip } from 'recharts';
import { PRESET_CELL_LINES } from '../constants';
import { Button } from './Button';
import { Input, RangeSlider } from './Input';
import { SweepLineChart } from './ResultsChart';
import { ASCENDService } from '../services/apiLive';
import {
  BatchSnapshotResponse,
  BatchSweepRequest,
  EnrichmentTerm,
  KeyGeneScore,
  MetabolomicsTaskStatusResponse,
  StringNetworkResponse,
} from '../types';

interface BatchWorkflowPageProps {
  isDark: boolean;
}

const delay = (ms: number) => new Promise((resolve) => setTimeout(resolve, ms));

const downloadBlob = (blob: Blob, filename: string) => {
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
};

export const BatchWorkflowPage: React.FC<BatchWorkflowPageProps> = ({ isDark }) => {
  const [smiles, setSmiles] = useState<string>('CC(=O)OC1=CC=CC=C1C(=O)O');
  const [cellLineId, setCellLineId] = useState<string>(PRESET_CELL_LINES[0].id);
  const [sweepVariable, setSweepVariable] = useState<'time' | 'dose'>('dose');
  const [rangeStart, setRangeStart] = useState<number>(0);
  const [rangeEnd, setRangeEnd] = useState<number>(100);
  const [rangeSteps, setRangeSteps] = useState<number>(6);
  const [fixedParamValue, setFixedParamValue] = useState<number>(24);
  const [snapshotValue, setSnapshotValue] = useState<number>(50);

  const [epochs, setEpochs] = useState<number>(100);
  const [imputation, setImputation] = useState<boolean>(false);
  const [nClusters, setNClusters] = useState<number>(4);
  const [enrichmentLibrary, setEnrichmentLibrary] = useState<string>('GO_Biological_Process_2021');
  const [enrichmentTopGenes, setEnrichmentTopGenes] = useState<number>(200);

  const [sweepResults, setSweepResults] = useState<any[] | null>(null);
  const [sweepTopGenes, setSweepTopGenes] = useState<string[]>([]);
  const [metabolomicsStatus, setMetabolomicsStatus] = useState<MetabolomicsTaskStatusResponse | null>(null);
  const [metabolomicsJobId, setMetabolomicsJobId] = useState<string | null>(null);
  const [snapshot, setSnapshot] = useState<BatchSnapshotResponse | null>(null);
  const [snapshotFlux, setSnapshotFlux] = useState<Record<string, number>>({});
  const [enrichment, setEnrichment] = useState<EnrichmentTerm[]>([]);
  const [ppi, setPpi] = useState<StringNetworkResponse | null>(null);
  const [keyGenes, setKeyGenes] = useState<KeyGeneScore[]>([]);

  const [loadingStep, setLoadingStep] = useState<'step1' | 'step2' | 'step3' | 'step4' | 'step5' | 'step6' | null>(null);
  const [errorMsg, setErrorMsg] = useState<string | null>(null);

  const sweepPayload: BatchSweepRequest = {
    smiles,
    cellLineId,
    sweepVariable,
    fixedParamValue,
    range: {
      start: rangeStart,
      end: rangeEnd,
      steps: rangeSteps,
    },
  };

  const enrichmentChart = useMemo(() => {
    return enrichment
      .slice(0, 12)
      .map((term) => ({
        term: term.term,
        score: -Math.log10((term.adjP ?? term.pValue) + 1e-12),
      }))
      .reverse();
  }, [enrichment]);

  const axisColor = isDark ? '#94a3b8' : '#64748b';
  const gridColor = isDark ? '#1e293b' : '#e2e8f0';

  const runStep1Sweep = async () => {
    setLoadingStep('step1');
    setErrorMsg(null);
    setMetabolomicsStatus(null);
    setMetabolomicsJobId(null);
    setSnapshot(null);
    setSnapshotFlux({});
    setEnrichment([]);
    setPpi(null);
    setKeyGenes([]);
    try {
      const response = await ASCENDService.predictBatchSweep(sweepPayload);
      if (response.type !== 'sweep' || !Array.isArray(response.data)) {
        throw new Error('Unexpected sweep response format.');
      }
      setSweepResults(response.data);
      setSweepTopGenes(response.topGenes ?? []);
      setSnapshotValue(rangeStart);
    } catch (err: any) {
      setSweepResults(null);
      setErrorMsg(err?.message ?? 'Step 1 failed.');
    } finally {
      setLoadingStep(null);
    }
  };

  const downloadAllTranscriptomics = async () => {
    try {
      const blob = await ASCENDService.fetchSweepMetabolomicsCsv(sweepPayload);
      const filename = `batch_transcriptomics_${sweepVariable}_${rangeStart}-${rangeEnd}_${rangeSteps}.csv`;
      downloadBlob(blob, filename);
    } catch (err: any) {
      setErrorMsg(err?.message ?? 'Failed to download transcriptomics matrix.');
    }
  };

  const runStep2Metabolomics = async () => {
    if (!sweepResults?.length) {
      setErrorMsg('Run Step 1 first.');
      return;
    }
    setLoadingStep('step2');
    setErrorMsg(null);
    try {
      const blob = await ASCENDService.fetchSweepMetabolomicsCsv(sweepPayload);
      const file = new File(
        [blob],
        `metabolomics_input_${sweepVariable}_${rangeStart}-${rangeEnd}_${rangeSteps}.csv`,
        { type: 'text/csv' },
      );
      const started = await ASCENDService.startMetabolomicsAnalysis({
        file,
        epochs,
        imputation,
        nClusters,
      });
      setMetabolomicsJobId(started.job_id);

      for (let i = 0; i < 180; i++) {
        const status = await ASCENDService.getMetabolomicsStatus(started.task_id);
        setMetabolomicsStatus(status);
        if (status.status === 'SUCCESS' || status.status === 'FAILURE') {
          break;
        }
        await delay(2000);
      }
    } catch (err: any) {
      setErrorMsg(err?.message ?? 'Step 2 failed.');
    } finally {
      setLoadingStep(null);
    }
  };

  const runStep3Snapshot = async () => {
    if (!sweepResults?.length) {
      setErrorMsg('Run Step 1 first.');
      return;
    }
    setLoadingStep('step3');
    setErrorMsg(null);
    setEnrichment([]);
    setPpi(null);
    setKeyGenes([]);
    try {
      const response = await ASCENDService.getBatchSnapshot({
        ...sweepPayload,
        snapshotValue,
      });
      setSnapshot(response);
      if (metabolomicsJobId) {
        const fluxResp = await ASCENDService.getBatchMetabolomicsSnapshot(metabolomicsJobId, response.pointIndex);
        setSnapshotFlux(fluxResp.flux);
      } else {
        setSnapshotFlux({});
      }
    } catch (err: any) {
      setSnapshot(null);
      setSnapshotFlux({});
      setErrorMsg(err?.message ?? 'Step 3 failed.');
    } finally {
      setLoadingStep(null);
    }
  };

  const runStep4Enrichment = async () => {
    if (!snapshot?.expression?.length) {
      setErrorMsg('Run Step 3 first.');
      return;
    }
    setLoadingStep('step4');
    setErrorMsg(null);
    try {
      const ranked = [...snapshot.expression].sort((a, b) => Math.abs(b.expressionLevel) - Math.abs(a.expressionLevel));
      const geneIds = ranked.slice(0, enrichmentTopGenes).map((g) => g.geneId);
      const response = await ASCENDService.runEnrichment({
        geneIds,
        library: enrichmentLibrary,
        topN: 20,
      });
      setEnrichment(response.data);
    } catch (err: any) {
      setEnrichment([]);
      setErrorMsg(err?.message ?? 'Step 4 failed.');
    } finally {
      setLoadingStep(null);
    }
  };

  const runStep5Ppi = async () => {
    if (!snapshot?.expression?.length) {
      setErrorMsg('Run Step 3 first.');
      return;
    }
    setLoadingStep('step5');
    setErrorMsg(null);
    try {
      const response = await ASCENDService.fetchStringNetwork({
        genes: snapshot.expression,
        minExpression: 1.5,
        species: 9606,
        requiredScore: 400,
        networkType: 'functional',
        maxGenes: 200,
      });
      setPpi(response);
    } catch (err: any) {
      setPpi(null);
      setErrorMsg(err?.message ?? 'Step 5 failed.');
    } finally {
      setLoadingStep(null);
    }
  };

  const runStep6KeyGenes = async () => {
    if (!snapshot?.expression?.length) {
      setErrorMsg('Run Step 3 first.');
      return;
    }
    setLoadingStep('step6');
    setErrorMsg(null);
    try {
      let flux = snapshotFlux;
      if ((!flux || Object.keys(flux).length === 0) && metabolomicsJobId) {
        const fluxResp = await ASCENDService.getBatchMetabolomicsSnapshot(metabolomicsJobId, snapshot.pointIndex);
        flux = fluxResp.flux;
        setSnapshotFlux(fluxResp.flux);
      }
      const response = await ASCENDService.identifyKeyGenes({
        expression: snapshot.expression,
        enrichment,
        ppiGenes: ppi?.genesUsed ?? [],
        metabolomicsFlux: flux ?? {},
        topN: 5,
      });
      setKeyGenes(response.topGenes);
    } catch (err: any) {
      setKeyGenes([]);
      setErrorMsg(err?.message ?? 'Step 6 failed.');
    } finally {
      setLoadingStep(null);
    }
  };

  const metabolomicsResult = metabolomicsStatus?.result;
  const imageKeys = [
    { key: 'umap_cluster', title: 'UMAP by Cluster' },
    { key: 'umap_stress', title: 'UMAP by Stress' },
    { key: 'heatmap_module', title: 'Cluster x Module' },
    { key: 'heatmap_pathway', title: 'Cluster x Pathway' },
    { key: 'volcano', title: 'Volcano Plot' },
    { key: 'stress_corr', title: 'Flux-Stress Correlation' },
  ] as const;

  const buildResultFileUrl = (fullPath: string) => {
    if (!metabolomicsJobId) return '#';
    const parts = fullPath.split(/[/\\]/);
    const filename = parts[parts.length - 1];
    return ASCENDService.getMetabolomicsResultUrl(metabolomicsJobId, filename);
  };

  return (
    <div className="space-y-8">
      <div className="text-center max-w-4xl mx-auto">
        <h1 className="text-4xl md:text-5xl font-bold text-slate-900 dark:text-white mb-3 tracking-tight">
          Batch <span className="text-transparent bg-clip-text bg-gradient-to-r from-cyan-600 to-blue-600 dark:from-cyan-400 dark:to-blue-500">Prediction Workflow</span>
        </h1>
        <p className="text-sm md:text-base text-slate-500 dark:text-slate-400">
          Dynamic time/dose sweep, metabolomics analysis, snapshot-driven enrichment and PPI, then integrated key-gene ranking.
        </p>
      </div>

      {errorMsg && (
        <div className="max-w-6xl mx-auto p-3 rounded-lg bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-500/40 text-sm text-red-600 dark:text-red-300">
          {errorMsg}
        </div>
      )}

      <div className="max-w-6xl mx-auto space-y-6">
        <section className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md">
          <div className="flex items-center justify-between gap-4 mb-4">
            <h2 className="text-lg font-bold text-slate-900 dark:text-white flex items-center gap-2">
              <Activity className="w-5 h-5 text-cyan-600 dark:text-cyan-400" /> Step 1: ASCEND Dynamic Transcriptomics Sweep
            </h2>
            <div className="flex gap-2">
              <Button onClick={runStep1Sweep} isLoading={loadingStep === 'step1'}>Run Step 1</Button>
              <Button variant="secondary" onClick={downloadAllTranscriptomics}>
                <Download className="w-4 h-4" /> Download All Genes
              </Button>
            </div>
          </div>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <Input label="Drug SMILES" value={smiles} onChange={(e) => setSmiles(e.target.value)} />
            <div>
              <label className="block mb-2 text-xs font-mono tracking-widest text-slate-500 dark:text-cyan-500 uppercase font-bold">Cell Line</label>
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
            <div>
              <label className="block mb-2 text-xs font-mono tracking-widest text-slate-500 dark:text-cyan-500 uppercase font-bold">Sweep Variable</label>
              <select
                value={sweepVariable}
                onChange={(e) => setSweepVariable(e.target.value as 'time' | 'dose')}
                className="w-full px-4 py-3 text-sm text-slate-900 dark:text-white bg-slate-50 dark:bg-slate-950/50 border border-slate-200 dark:border-slate-700 rounded-lg"
              >
                <option value="dose">Dose Sweep</option>
                <option value="time">Time Sweep</option>
              </select>
            </div>
          </div>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mt-2">
            <Input label={sweepVariable === 'dose' ? 'Fixed Time (h)' : 'Fixed Dose (uM)'} type="number" value={fixedParamValue} onChange={(e) => setFixedParamValue(Number(e.target.value))} />
            <Input label="Range Start" type="number" value={rangeStart} onChange={(e) => setRangeStart(Number(e.target.value))} />
            <Input label="Range End" type="number" value={rangeEnd} onChange={(e) => setRangeEnd(Number(e.target.value))} />
          </div>
          <Input label="Range Steps" type="number" min={2} value={rangeSteps} onChange={(e) => setRangeSteps(Number(e.target.value))} />
          {sweepTopGenes.length > 0 ? (
            <p className="text-xs text-slate-500 dark:text-slate-400 mb-3">Top variable genes used in curves: {sweepTopGenes.join(', ')}</p>
          ) : null}
          {sweepResults?.length ? (
            <SweepLineChart
              data={sweepResults as any}
              isDark={isDark}
              xAxisLabel={sweepVariable === 'dose' ? 'Dose (uM)' : 'Time (h)'}
            />
          ) : null}
        </section>

        <section className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md">
          <div className="flex items-center justify-between gap-4 mb-4">
            <h2 className="text-lg font-bold text-slate-900 dark:text-white flex items-center gap-2">
              <Beaker className="w-5 h-5 text-cyan-600 dark:text-cyan-400" /> Step 2: Metabolomics Prediction and Plotting
            </h2>
            <Button onClick={runStep2Metabolomics} isLoading={loadingStep === 'step2'}>Run Step 2</Button>
          </div>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <Input label="Epochs" type="number" min={10} value={epochs} onChange={(e) => setEpochs(Number(e.target.value))} />
            <Input label="Clusters" type="number" min={2} value={nClusters} onChange={(e) => setNClusters(Number(e.target.value))} />
            <div className="flex items-end pb-4">
              <label className="inline-flex items-center gap-2 text-sm text-slate-700 dark:text-slate-300">
                <input type="checkbox" checked={imputation} onChange={(e) => setImputation(e.target.checked)} />
                Enable Imputation
              </label>
            </div>
          </div>
          {metabolomicsStatus ? (
            <div className="space-y-4">
              <div className="text-sm text-slate-600 dark:text-slate-300">
                Status: <span className="font-semibold">{metabolomicsStatus.status}</span> ({metabolomicsStatus.progress}%)
                <div>{metabolomicsStatus.message}</div>
              </div>
              {metabolomicsResult?.files ? (
                <div className="flex flex-wrap gap-3">
                  {Object.entries(metabolomicsResult.files).map(([label, filePath]) => (
                    <a key={label} href={buildResultFileUrl(filePath)} target="_blank" rel="noreferrer">
                      <Button variant="secondary" className="flex items-center gap-2">
                        <Download className="w-4 h-4" /> Download {label}
                      </Button>
                    </a>
                  ))}
                </div>
              ) : null}
              {metabolomicsResult?.images ? (
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                  {imageKeys.map((item) => {
                    const filename = (metabolomicsResult.images as any)[item.key];
                    if (!filename || !metabolomicsJobId) return null;
                    return (
                      <div key={item.key} className="border border-slate-200 dark:border-slate-800 rounded-lg p-3">
                        <p className="text-sm font-semibold text-slate-900 dark:text-white mb-2">{item.title}</p>
                        <img src={ASCENDService.getMetabolomicsResultUrl(metabolomicsJobId, filename)} alt={item.title} className="w-full h-auto rounded" />
                      </div>
                    );
                  })}
                </div>
              ) : null}
            </div>
          ) : null}
        </section>

        <section className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md">
          <div className="flex items-center justify-between gap-4 mb-4">
            <h2 className="text-lg font-bold text-slate-900 dark:text-white flex items-center gap-2">
              <FlaskConical className="w-5 h-5 text-cyan-600 dark:text-cyan-400" /> Step 3: Select Snapshot Point
            </h2>
            <Button onClick={runStep3Snapshot} isLoading={loadingStep === 'step3'}>Run Step 3</Button>
          </div>
          <RangeSlider
            label={sweepVariable === 'dose' ? 'Snapshot Dose' : 'Snapshot Time'}
            min={Math.min(rangeStart, rangeEnd)}
            max={Math.max(rangeStart, rangeEnd)}
            unit={sweepVariable === 'dose' ? 'uM' : 'h'}
            value={snapshotValue}
            onChange={(e) => setSnapshotValue(Number(e.target.value))}
          />
          {snapshot ? (
            <p className="text-sm text-slate-600 dark:text-slate-300">
              Snapshot selected: x={snapshot.xValue}, time={snapshot.time}h, dose={snapshot.dose}uM, point={snapshot.cellId}
            </p>
          ) : null}
        </section>

        <section className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md">
          <div className="flex items-center justify-between gap-4 mb-4">
            <h2 className="text-lg font-bold text-slate-900 dark:text-white">Step 4: Enrichment Analysis</h2>
            <Button onClick={runStep4Enrichment} isLoading={loadingStep === 'step4'}>Run Step 4</Button>
          </div>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-4">
            <Input label="Top Differential Genes" type="number" min={20} max={1000} value={enrichmentTopGenes} onChange={(e) => setEnrichmentTopGenes(Number(e.target.value))} />
            <div>
              <label className="block mb-2 text-xs font-mono tracking-widest text-slate-500 dark:text-cyan-500 uppercase font-bold">GO Library</label>
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
              <Network className="w-5 h-5 text-cyan-600 dark:text-cyan-400" /> Step 5: PPI Analysis
            </h2>
            <Button onClick={runStep5Ppi} isLoading={loadingStep === 'step5'}>Run Step 5</Button>
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
              <KeyRound className="w-5 h-5 text-cyan-600 dark:text-cyan-400" /> Step 6: Identify 5 Key Genes
            </h2>
            <Button onClick={runStep6KeyGenes} isLoading={loadingStep === 'step6'}>Run Step 6</Button>
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
