import React, { useEffect, useMemo, useState } from 'react';
import {
  Activity,
  FlaskConical,
  Image,
  Layers,
  Upload,
  Download
} from 'lucide-react';
import {
  ResponsiveContainer,
  BarChart,
  Bar,
  XAxis,
  YAxis,
  Tooltip,
  CartesianGrid
} from 'recharts';
import { Button } from './Button';
import { Input } from './Input';
import { ASCENDService } from '../services/apiLive';
import {
  MetabolomicsTaskStatusResponse,
  MetabolomicsTaskResponse,
  MetabolomicsDiffPlotResponse
} from '../types';

interface MetabolomicsPageProps {
  isDark: boolean;
  transcriptomicsSeedFile?: File | null;
  transcriptomicsSeedLabel?: string | null;
}

export const MetabolomicsPage: React.FC<MetabolomicsPageProps> = ({
  isDark,
  transcriptomicsSeedFile,
  transcriptomicsSeedLabel
}) => {
  const [file, setFile] = useState<File | null>(null);
  const [epochs, setEpochs] = useState<number>(100);
  const [imputation, setImputation] = useState<boolean>(false);
  const [nClusters, setNClusters] = useState<number>(4);

  const [jobId, setJobId] = useState<string | null>(null);
  const [taskId, setTaskId] = useState<string | null>(null);
  const [status, setStatus] = useState<MetabolomicsTaskStatusResponse | null>(null);

  const [isStarting, setIsStarting] = useState(false);
  const [errorMsg, setErrorMsg] = useState<string | null>(null);

  const [diffC1, setDiffC1] = useState<number>(0);
  const [diffC2, setDiffC2] = useState<number>(1);
  const [diffLoading, setDiffLoading] = useState<boolean>(false);
  const [customVolcano, setCustomVolcano] = useState<string | null>(null);
  const [customCsv, setCustomCsv] = useState<string | null>(null);
  const [prefillLabel, setPrefillLabel] = useState<string | null>(null);
  const [fluxRows, setFluxRows] = useState<Record<string, Record<string, number>> | null>(null);
  const [fluxCells, setFluxCells] = useState<string[]>([]);
  const [selectedCell, setSelectedCell] = useState<string>('');
  const [fluxLoading, setFluxLoading] = useState<boolean>(false);
  const [fluxError, setFluxError] = useState<string | null>(null);

  const handleFileChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const selected = e.target.files?.[0] ?? null;
    setFile(selected);
    setErrorMsg(null);
  };

  useEffect(() => {
    if (transcriptomicsSeedFile && !file) {
      setFile(transcriptomicsSeedFile);
      setPrefillLabel(transcriptomicsSeedLabel ?? transcriptomicsSeedFile.name);
      setErrorMsg(null);
    }
  }, [transcriptomicsSeedFile, transcriptomicsSeedLabel, file]);

  useEffect(() => {
    const loadFlux = async () => {
      if (status?.status !== 'SUCCESS' || !jobId) {
        return;
      }
      setFluxLoading(true);
      setFluxError(null);
      try {
        const response = await fetch(getResultUrl('predicted_flux.csv'));
        if (!response.ok) {
          throw new Error(`Failed to load flux data (${response.status})`);
        }
        const csvText = await response.text();
        const lines = csvText.trim().split(/\r?\n/).filter(Boolean);
        if (lines.length < 2) {
          throw new Error('Flux data is empty');
        }
        const headers = lines[0].split(',');
        const modules = headers.slice(1).map((h) => h.replace(/^\"|\"$/g, ''));
        const rows: Record<string, Record<string, number>> = {};
        const cells: string[] = [];
        for (let i = 1; i < lines.length; i++) {
          const parts = lines[i].split(',');
          if (parts.length < 2) continue;
          const cell = parts[0] || `cell_${i}`;
          const row: Record<string, number> = {};
          for (let j = 1; j < parts.length && j <= modules.length; j++) {
            const value = parseFloat(parts[j]);
            row[modules[j - 1]] = Number.isFinite(value) ? value : 0;
          }
          rows[cell] = row;
          cells.push(cell);
        }
        setFluxRows(rows);
        setFluxCells(cells);
        setSelectedCell((prev) => prev || cells[0] || '');
      } catch (err: any) {
        setFluxError(err?.message ?? 'Failed to load flux data');
      } finally {
        setFluxLoading(false);
      }
    };
    loadFlux();
  }, [status?.status, jobId]);

  const handleStart = async () => {
    if (!file) {
      setErrorMsg('Please upload a transcriptome CSV file first.');
      return;
    }
    setIsStarting(true);
    setErrorMsg(null);
    setStatus(null);
    setCustomVolcano(null);
    setCustomCsv(null);
    try {
      const response: MetabolomicsTaskResponse = await ASCENDService.startMetabolomicsAnalysis({
        file,
        epochs,
        imputation,
        nClusters
      });
      setJobId(response.job_id);
      setTaskId(response.task_id);
    } catch (err: any) {
      setErrorMsg(err?.message ?? 'Failed to start metabolomics analysis.');
    } finally {
      setIsStarting(false);
    }
  };

  const handleDiffPlot = async () => {
    if (!jobId) {
      setErrorMsg('No job available. Run an analysis first.');
      return;
    }
    setDiffLoading(true);
    setErrorMsg(null);
    try {
      const response: MetabolomicsDiffPlotResponse = await ASCENDService.generateMetabolomicsDiffPlot(
        jobId,
        diffC1,
        diffC2
      );
      if (response.status === 'success') {
        setCustomVolcano(response.image);
        setCustomCsv(response.csv);
      } else {
        setErrorMsg(response.error ?? 'Failed to generate differential plot.');
      }
    } catch (err: any) {
      setErrorMsg(err?.message ?? 'Failed to generate differential plot.');
    } finally {
      setDiffLoading(false);
    }
  };

  useEffect(() => {
    if (!taskId) {
      return;
    }
    let isActive = true;
    const poll = async () => {
      try {
        const nextStatus = await ASCENDService.getMetabolomicsStatus(taskId);
        if (!isActive) return;
        setStatus(nextStatus);
        if (['SUCCESS', 'FAILURE', 'UNKNOWN'].includes(nextStatus.status)) {
          isActive = false;
        }
      } catch (err: any) {
        if (!isActive) return;
        setErrorMsg(err?.message ?? 'Failed to retrieve task status.');
      }
    };
    const interval = setInterval(poll, 2000);
    poll();
    return () => {
      isActive = false;
      clearInterval(interval);
    };
  }, [taskId]);

  const progress = status?.progress ?? 0;
  const statusLabel = status?.status ?? 'IDLE';
  const resultImages = status?.result?.images;

  const defaultVolcano = customVolcano ?? resultImages?.volcano;
  const volcanoTitle = customVolcano ? `Custom Volcano (C${diffC1} vs C${diffC2})` : 'Default Volcano (C0 vs C1)';

  const getResultUrl = (filename?: string) => {
    if (!filename || !jobId) return '';
    return ASCENDService.getMetabolomicsResultUrl(jobId, filename);
  };

  const imageCards = useMemo(() => {
    if (!resultImages) return [];
    return [
      { key: 'umap_cluster', title: 'UMAP by Cluster', description: 'Metabolic flux embedding colored by cluster.' },
      { key: 'umap_stress', title: 'UMAP by Stress', description: 'Stress score mapped onto UMAP space.' },
      { key: 'heatmap_module', title: 'Cluster x Module', description: 'Cluster-wise module activity heatmap.' },
      { key: 'heatmap_pathway', title: 'Cluster x Pathway', description: 'Aggregated pathway activity heatmap.' },
      { key: 'stress_corr', title: 'Flux-Stress Correlation', description: 'Top metabolic modules correlated with stress.' }
    ] as const;
  }, [resultImages]);

  const fluxChartData = useMemo(() => {
    if (!fluxRows || !selectedCell) return [];
    const row = fluxRows[selectedCell];
    if (!row) return [];
    const entries = Object.entries(row).map(([module, value]) => ({
      module,
      value
    }));
    entries.sort((a, b) => Math.abs(b.value) - Math.abs(a.value));
    return entries.slice(0, 20).reverse();
  }, [fluxRows, selectedCell]);

  const axisColor = isDark ? '#94a3b8' : '#64748b';
  const gridColor = isDark ? '#1e293b' : '#e2e8f0';

  return (
    <>
      <div className="mb-10 text-center max-w-3xl mx-auto">
        <h1 className="text-4xl md:text-5xl font-bold text-slate-900 dark:text-white mb-4 tracking-tight transition-colors">
          Predict <span className="text-transparent bg-clip-text bg-gradient-to-r from-cyan-600 to-blue-600 dark:from-cyan-400 dark:to-blue-500">Metabolic Flux</span>
        </h1>
        <p className="text-sm md:text-base text-slate-500 dark:text-slate-400">
          Upload transcriptome data, train the scFEA model, and explore pathway-level metabolomics predictions.
        </p>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-12 gap-8">
        <div className="lg:col-span-4 space-y-6">
          <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-xl shadow-slate-200/40 dark:shadow-none relative overflow-hidden transition-all duration-300">
            <div className="absolute top-0 left-0 w-full h-1 bg-gradient-to-r from-cyan-500 to-blue-600"></div>
            <h2 className="text-lg font-bold text-slate-900 dark:text-white mb-6 flex items-center gap-2">
              <Layers className="w-5 h-5 text-cyan-600 dark:text-cyan-400" />
              Transcriptome Upload
            </h2>

            <div className="border-2 border-dashed border-slate-200 dark:border-slate-700 rounded-xl p-6 text-center bg-slate-50 dark:bg-slate-950/50">
              <input
                id="metabo-file"
                type="file"
                accept=".csv"
                onChange={handleFileChange}
                className="hidden"
              />
              <label htmlFor="metabo-file" className="flex flex-col items-center gap-3 cursor-pointer">
                <div className="p-3 rounded-full bg-cyan-50 dark:bg-cyan-950/40 border border-cyan-200 dark:border-cyan-700">
                  <Upload className="w-5 h-5 text-cyan-600 dark:text-cyan-400" />
                </div>
                <div>
                  <p className="text-sm font-semibold text-slate-800 dark:text-slate-200">
                    {file ? file.name : 'Select transcriptome CSV'}
                  </p>
                  <p className="text-xs text-slate-500 dark:text-slate-400">
                    {file ? `${(file.size / 1024).toFixed(1)} KB` : 'CSV with genes as columns, cells as rows'}
                  </p>
                </div>
              </label>
            </div>
            {transcriptomicsSeedFile && (
              <div className="mt-4 flex flex-col sm:flex-row sm:items-center sm:justify-between gap-3 rounded-lg border border-cyan-200/50 dark:border-cyan-700/40 bg-cyan-50/60 dark:bg-cyan-950/20 p-3 text-sm text-slate-700 dark:text-slate-300">
                <div>
                  <p className="font-semibold text-cyan-700 dark:text-cyan-300">Transcriptomics output available</p>
                  <p className="text-xs text-slate-500 dark:text-slate-400">{prefillLabel ?? transcriptomicsSeedFile.name}</p>
                </div>
                <Button
                  variant="outline"
                  className="flex items-center gap-2"
                  onClick={() => {
                    setFile(transcriptomicsSeedFile);
                    setPrefillLabel(transcriptomicsSeedLabel ?? transcriptomicsSeedFile.name);
                  }}
                >
                  <Layers className="w-4 h-4" /> Use Output
                </Button>
              </div>
            )}

            <div className="mt-6">
              <Input
                label="Training Epochs"
                type="number"
                min={10}
                max={500}
                value={epochs}
                onChange={(e) => {
                  const value = Number(e.target.value);
                  if (!Number.isNaN(value)) setEpochs(value);
                }}
              />
              <Input
                label="Cluster Count"
                type="number"
                min={2}
                max={12}
                value={nClusters}
                onChange={(e) => {
                  const value = Number(e.target.value);
                  if (!Number.isNaN(value)) setNClusters(value);
                }}
              />
              <label className="flex items-center gap-3 text-sm text-slate-600 dark:text-slate-300 mt-2">
                <input
                  type="checkbox"
                  checked={imputation}
                  onChange={(e) => setImputation(e.target.checked)}
                  className="h-4 w-4 accent-cyan-500"
                />
                Enable MAGIC imputation for sparse single-cell profiles
              </label>
            </div>

            <div className="mt-6">
              <Button onClick={handleStart} isLoading={isStarting} className="w-full justify-center">
                <FlaskConical className="w-4 h-4" />
                Start Metabolomics Prediction
              </Button>
            </div>

            {errorMsg && (
              <div className="mt-4 p-3 rounded-lg bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-500/40 text-sm text-red-600 dark:text-red-300">
                {errorMsg}
              </div>
            )}
          </div>

          <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md">
            <h3 className="text-sm font-bold text-slate-900 dark:text-white mb-4 flex items-center gap-2">
              <Activity className="w-4 h-4 text-cyan-600 dark:text-cyan-400" />
              Task Status
            </h3>
            <div className="flex items-center justify-between text-xs text-slate-500 dark:text-slate-400 font-mono">
              <span>Status</span>
              <span>{statusLabel}</span>
            </div>
            <div className="mt-3 h-2 w-full bg-slate-200 dark:bg-slate-800 rounded-full overflow-hidden">
              <div
                className="h-full bg-gradient-to-r from-cyan-500 to-blue-500 transition-all duration-500"
                style={{ width: `${Math.min(progress, 100)}%` }}
              ></div>
            </div>
            <p className="mt-3 text-xs text-slate-500 dark:text-slate-400">
              {status?.message ?? 'Upload a dataset to begin.'}
            </p>
            {jobId && (
              <p className="mt-2 text-[11px] text-slate-400 dark:text-slate-500 font-mono break-all">
                Job ID: {jobId}
              </p>
            )}
          </div>
        </div>

        <div className="lg:col-span-8 space-y-6">
          <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md">
            <h2 className="text-lg font-bold text-slate-900 dark:text-white mb-3 flex items-center gap-2">
              <Image className="w-5 h-5 text-cyan-600 dark:text-cyan-400" />
              Metabolomics Outputs
            </h2>
            <p className="text-sm text-slate-500 dark:text-slate-400">
              Visual summaries appear here once training and downstream analysis complete.
            </p>
          </div>

          {status?.status === 'SUCCESS' && (
            <>
              <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md">
                <div className="flex flex-col md:flex-row md:items-center md:justify-between gap-4 mb-4">
                  <div>
                    <h3 className="text-lg font-bold text-slate-900 dark:text-white">Numerical Flux Summary</h3>
                    <p className="text-sm text-slate-500 dark:text-slate-400">
                      Top metabolic modules for the selected cell. Values are predicted flux scores.
                    </p>
                  </div>
                  {fluxCells.length > 1 && (
                    <div className="min-w-[200px]">
                      <label className="block mb-2 text-xs font-mono tracking-widest text-slate-500 dark:text-cyan-500 uppercase font-bold">
                        Cell
                      </label>
                      <select
                        value={selectedCell}
                        onChange={(e) => setSelectedCell(e.target.value)}
                        className="w-full appearance-none bg-slate-50 dark:bg-slate-950/50 border border-slate-200 dark:border-slate-700 rounded-lg px-4 py-3 text-sm text-slate-900 dark:text-white outline-none"
                      >
                        {fluxCells.map((cell) => (
                          <option key={cell} value={cell}>{cell}</option>
                        ))}
                      </select>
                    </div>
                  )}
                </div>

                {fluxLoading && (
                  <p className="text-sm text-slate-500 dark:text-slate-400">Loading flux data...</p>
                )}
                {fluxError && (
                  <div className="p-3 rounded-lg bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-500/40 text-sm text-red-600 dark:text-red-300">
                    {fluxError}
                  </div>
                )}
                {!fluxLoading && !fluxError && fluxChartData.length > 0 && (
                  <div className="h-[320px]">
                    <ResponsiveContainer width="100%" height="100%">
                      <BarChart data={fluxChartData} layout="vertical" margin={{ left: 8, right: 24 }}>
                        <CartesianGrid stroke={gridColor} strokeDasharray="3 3" />
                        <XAxis type="number" stroke={axisColor} />
                        <YAxis dataKey="module" type="category" width={80} stroke={axisColor} />
                        <Tooltip />
                        <Bar dataKey="value" fill="#06b6d4" />
                      </BarChart>
                    </ResponsiveContainer>
                  </div>
                )}
                {!fluxLoading && !fluxError && fluxChartData.length === 0 && (
                  <p className="text-sm text-slate-500 dark:text-slate-400">
                    No flux data available for charting.
                  </p>
                )}
              </div>

              {resultImages && (
              <>
              <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md">
                <div className="flex flex-col md:flex-row md:items-center md:justify-between gap-4 mb-4">
                  <div>
                    <h3 className="text-lg font-bold text-slate-900 dark:text-white">Differential Flux Plot</h3>
                    <p className="text-sm text-slate-500 dark:text-slate-400">
                      Generate a volcano plot comparing any two clusters.
                    </p>
                  </div>
                  <div className="flex flex-wrap gap-3">
                    <div className="w-24">
                      <Input
                        label="Cluster A"
                        type="number"
                        min={0}
                        value={diffC1}
                        onChange={(e) => setDiffC1(Number(e.target.value))}
                      />
                    </div>
                    <div className="w-24">
                      <Input
                        label="Cluster B"
                        type="number"
                        min={0}
                        value={diffC2}
                        onChange={(e) => setDiffC2(Number(e.target.value))}
                      />
                    </div>
                    <div className="flex items-end">
                      <Button onClick={handleDiffPlot} isLoading={diffLoading} variant="outline">
                        Generate
                      </Button>
                    </div>
                  </div>
                </div>
                {defaultVolcano && (
                  <div className="border border-slate-200 dark:border-slate-800 rounded-xl overflow-hidden bg-slate-50 dark:bg-slate-950/50">
                    <div className="px-4 py-3 border-b border-slate-200 dark:border-slate-800 text-sm font-semibold text-slate-700 dark:text-slate-200">
                      {volcanoTitle}
                    </div>
                    <img
                      src={getResultUrl(defaultVolcano)}
                      alt={volcanoTitle}
                      className="w-full object-cover"
                    />
                  </div>
                )}
                <div className="mt-4 flex flex-wrap gap-3">
                  {customCsv && (
                    <a href={getResultUrl(customCsv)} download className="inline-flex">
                      <Button variant="secondary" className="flex items-center gap-2">
                        <Download className="w-4 h-4" /> Download Diff CSV
                      </Button>
                    </a>
                  )}
                  {defaultVolcano && (
                    <a href={getResultUrl(defaultVolcano)} download className="inline-flex">
                      <Button variant="secondary" className="flex items-center gap-2">
                        <Download className="w-4 h-4" /> Download Volcano PNG
                      </Button>
                    </a>
                  )}
                </div>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                {imageCards.map((card) => (
                  <div
                    key={card.key}
                    className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl overflow-hidden shadow-md"
                  >
                    <div className="px-4 py-3 border-b border-slate-200 dark:border-slate-800">
                      <h4 className="text-sm font-bold text-slate-800 dark:text-slate-100">{card.title}</h4>
                      <p className="text-xs text-slate-500 dark:text-slate-400">{card.description}</p>
                    </div>
                    <img
                      src={getResultUrl(resultImages[card.key])}
                      alt={card.title}
                      className="w-full object-cover"
                    />
                  </div>
                ))}
              </div>
              </>
              )}

              <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md">
                <h3 className="text-lg font-bold text-slate-900 dark:text-white mb-3">Download Data</h3>
                <div className="flex flex-wrap gap-3">
                  <a href={getResultUrl('predicted_flux.csv')} download className="inline-flex">
                    <Button variant="secondary" className="flex items-center gap-2">
                      <Download className="w-4 h-4" /> Predicted Flux
                    </Button>
                  </a>
                  <a href={getResultUrl('predicted_balance.csv')} download className="inline-flex">
                    <Button variant="secondary" className="flex items-center gap-2">
                      <Download className="w-4 h-4" /> Predicted Balance
                    </Button>
                  </a>
                  <a href={getResultUrl('PCA_coordinates.csv')} download className="inline-flex">
                    <Button variant="secondary" className="flex items-center gap-2">
                      <Download className="w-4 h-4" /> PCA Coordinates
                    </Button>
                  </a>
                  <a href={getResultUrl('cell_clusters.csv')} download className="inline-flex">
                    <Button variant="secondary" className="flex items-center gap-2">
                      <Download className="w-4 h-4" /> Cluster Labels
                    </Button>
                  </a>
                </div>
              </div>
            </>
          )}

          {status?.status === 'FAILURE' && (
            <div className="p-4 rounded-xl bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-500/40 text-sm text-red-600 dark:text-red-300">
              {status.error ?? 'Metabolomics analysis failed. Check the backend logs for details.'}
            </div>
          )}
        </div>
      </div>
    </>
  );
};
