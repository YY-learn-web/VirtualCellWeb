import React, { useMemo, useState } from 'react';
import {
  BarChart3,
  Database,
  Download,
  FlaskConical,
  Share2
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
import { ASCENDService } from '../services/api2';
import { GeneExpressionResult, EnrichmentResponse, EnrichmentTerm } from '../types';

interface EnrichmentPageProps {
  isDark: boolean;
  differentialResults: GeneExpressionResult[] | null;
  contextLabel?: string;
}

export const EnrichmentPage: React.FC<EnrichmentPageProps> = ({
  isDark,
  differentialResults,
  contextLabel
}) => {
  const [direction, setDirection] = useState<'up' | 'down' | 'both'>('both');
  const [topGenes, setTopGenes] = useState<number>(200);
  const [library, setLibrary] = useState<string>('GO_Biological_Process_2021');
  const [results, setResults] = useState<EnrichmentTerm[]>([]);
  const [loading, setLoading] = useState(false);
  const [errorMsg, setErrorMsg] = useState<string | null>(null);
  const [ppiImage, setPpiImage] = useState<string | null>(null);
  const [ppiLink, setPpiLink] = useState<string | null>(null);
  const [ppiStats, setPpiStats] = useState<{ filtered: number; mapped: number } | null>(null);
  const [ppiLoading, setPpiLoading] = useState(false);
  const [ppiError, setPpiError] = useState<string | null>(null);

  const availableLibraries = [
    'GO_Biological_Process_2021',
    'GO_Molecular_Function_2021',
    'GO_Cellular_Component_2021',
    'KEGG_2021_Human'
  ];

  const rankedGenes = useMemo(() => {
    if (!differentialResults?.length) return [];
    const filtered = differentialResults.filter((g) => {
      if (direction === 'up') return g.expressionLevel > 0;
      if (direction === 'down') return g.expressionLevel < 0;
      return true;
    });
    return [...filtered].sort((a, b) => Math.abs(b.expressionLevel) - Math.abs(a.expressionLevel));
  }, [differentialResults, direction]);

  const geneList = useMemo(() => {
    return rankedGenes.slice(0, topGenes).map((g) => g.geneId);
  }, [rankedGenes, topGenes]);

  const chartData = useMemo(() => {
    return results.map((term) => ({
      term: term.term,
      score: -Math.log10((term.adjP ?? term.pValue) + 1e-12)
    })).slice(0, 15).reverse();
  }, [results]);

  const axisColor = isDark ? '#94a3b8' : '#64748b';
  const gridColor = isDark ? '#1e293b' : '#e2e8f0';

  const runEnrichment = async () => {
    if (!geneList.length) {
      setErrorMsg('No differential genes available. Run transcriptomics prediction first.');
      return;
    }
    setLoading(true);
    setErrorMsg(null);
    try {
      const response: EnrichmentResponse = await ASCENDService.runEnrichment({
        geneIds: geneList,
        library,
        topN: 20
      });
      setResults(response.data);
    } catch (err: any) {
      setErrorMsg(err?.message ?? 'Enrichment failed.');
      setResults([]);
    } finally {
      setLoading(false);
    }
  };

  const runStringNetwork = async () => {
    if (!differentialResults?.length) {
      setPpiError('No transcriptomics results available. Run prediction first.');
      return;
    }
    setPpiLoading(true);
    setPpiError(null);
    try {
      const response = await ASCENDService.fetchStringNetwork({
        genes: differentialResults,
        minExpression: 1.5,
        species: 9606,
        requiredScore: 400,
        networkType: 'functional',
        maxGenes: 200
      });
      setPpiImage(response.image);
      setPpiLink(response.link ?? null);
      setPpiStats({ filtered: response.filteredCount, mapped: response.mappedCount });
    } catch (err: any) {
      setPpiError(err?.message ?? 'STRING network failed.');
      setPpiImage(null);
      setPpiLink(null);
      setPpiStats(null);
    } finally {
      setPpiLoading(false);
    }
  };

  const downloadResults = () => {
    if (!results.length) return;
    const header = ['term', 'pValue', 'adjP', 'combinedScore', 'genes'].join(',');
    const rows = results.map((r) => [
      `"${r.term.replace(/\"/g, '""')}"`,
      r.pValue,
      r.adjP ?? '',
      r.combinedScore ?? '',
      `"${(r.genes ?? []).join(';')}"`
    ].join(','));
    const csv = [header, ...rows].join('\n');
    const file = new File([csv], 'gene_enrichment_results.csv', { type: 'text/csv' });
    const url = URL.createObjectURL(file);
    const a = document.createElement('a');
    a.href = url;
    a.download = file.name;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  };

  return (
    <>
      <div className="mb-10 text-center max-w-3xl mx-auto">
        <h1 className="text-4xl md:text-5xl font-bold text-slate-900 dark:text-white mb-4 tracking-tight transition-colors">
          Gene <span className="text-transparent bg-clip-text bg-gradient-to-r from-cyan-600 to-blue-600 dark:from-cyan-400 dark:to-blue-500">Enrichment</span>
        </h1>
        <p className="text-sm md:text-base text-slate-500 dark:text-slate-400">
          Use differential transcriptomic predictions to uncover enriched biological pathways.
        </p>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-12 gap-8">
        <div className="lg:col-span-4 space-y-6">
          <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-xl shadow-slate-200/40 dark:shadow-none relative overflow-hidden transition-all duration-300">
            <div className="absolute top-0 left-0 w-full h-1 bg-gradient-to-r from-cyan-500 to-blue-600"></div>
            <h2 className="text-lg font-bold text-slate-900 dark:text-white mb-6 flex items-center gap-2">
              <Database className="w-5 h-5 text-cyan-600 dark:text-cyan-400" />
              Enrichment Setup
            </h2>

            <div className="mb-4">
              <label className="block mb-2 text-xs font-mono tracking-widest text-slate-500 dark:text-cyan-500 uppercase font-bold">
                Direction
              </label>
              <div className="flex gap-2">
                {(['both', 'up', 'down'] as const).map((opt) => (
                  <button
                    key={opt}
                    type="button"
                    onClick={() => setDirection(opt)}
                    className={`flex-1 py-2 text-xs rounded-lg transition-all border font-semibold ${
                      direction === opt
                        ? 'bg-cyan-50 dark:bg-cyan-900/30 text-cyan-700 dark:text-cyan-300 border-cyan-200 dark:border-cyan-700'
                        : 'bg-slate-50 dark:bg-slate-800/50 text-slate-500 dark:text-slate-400 border-slate-200 dark:border-slate-700'
                    }`}
                  >
                    {opt.toUpperCase()}
                  </button>
                ))}
              </div>
            </div>

            <Input
              label="Top Genes"
              type="number"
              min={20}
              max={1000}
              value={topGenes}
              onChange={(e) => setTopGenes(Number(e.target.value))}
            />

            <label className="block mb-2 text-xs font-mono tracking-widest text-slate-500 dark:text-cyan-500 uppercase font-bold">
              Gene Set Library
            </label>
            <select
              value={library}
              onChange={(e) => setLibrary(e.target.value)}
              className="w-full appearance-none bg-slate-50 dark:bg-slate-950/50 border border-slate-200 dark:border-slate-700 rounded-lg px-4 py-3 text-sm text-slate-900 dark:text-white outline-none"
            >
              {availableLibraries.map((lib) => (
                <option key={lib} value={lib}>{lib}</option>
              ))}
            </select>

            <div className="mt-6">
              <Button onClick={runEnrichment} isLoading={loading} className="w-full justify-center">
                <FlaskConical className="w-4 h-4" /> Run Enrichment
              </Button>
            </div>

            {contextLabel && (
              <p className="mt-4 text-xs text-slate-500 dark:text-slate-400">
                Source: {contextLabel}
              </p>
            )}

            {errorMsg && (
              <div className="mt-4 p-3 rounded-lg bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-500/40 text-sm text-red-600 dark:text-red-300">
                {errorMsg}
              </div>
            )}
          </div>
        </div>

        <div className="lg:col-span-8 space-y-6">
          <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md">
            <h2 className="text-lg font-bold text-slate-900 dark:text-white mb-3 flex items-center gap-2">
              <BarChart3 className="w-5 h-5 text-cyan-600 dark:text-cyan-400" />
              Enrichment Results
            </h2>
            <p className="text-sm text-slate-500 dark:text-slate-400">
              Top enriched terms ranked by adjusted p-value. Bars show -log10(adj P).
            </p>
          </div>

          {results.length > 0 && (
            <>
              <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md">
                <div className="flex items-center justify-between mb-4">
                  <h3 className="text-sm font-bold text-slate-900 dark:text-white">Top Enriched Terms</h3>
                  <Button variant="secondary" className="flex items-center gap-2" onClick={downloadResults}>
                    <Download className="w-4 h-4" /> Download CSV
                  </Button>
                </div>
                <div className="h-[320px]">
                  <ResponsiveContainer width="100%" height="100%">
                    <BarChart data={chartData} layout="vertical" margin={{ left: 12, right: 24 }}>
                      <CartesianGrid stroke={gridColor} strokeDasharray="3 3" />
                      <XAxis type="number" stroke={axisColor} />
                      <YAxis dataKey="term" type="category" width={140} stroke={axisColor} />
                      <Tooltip />
                      <Bar dataKey="score" fill="#06b6d4" />
                    </BarChart>
                  </ResponsiveContainer>
                </div>
              </div>

              <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md overflow-x-auto">
                <table className="min-w-full text-left text-sm">
                  <thead className="bg-slate-50 dark:bg-slate-800 text-slate-500 dark:text-slate-400 uppercase text-xs">
                    <tr>
                      <th className="px-4 py-2">Term</th>
                      <th className="px-4 py-2">Adj P</th>
                      <th className="px-4 py-2">Combined</th>
                      <th className="px-4 py-2">Genes</th>
                    </tr>
                  </thead>
                  <tbody className="text-slate-700 dark:text-slate-300">
                    {results.map((term) => (
                      <tr key={term.term} className="border-t border-slate-200 dark:border-slate-800">
                        <td className="px-4 py-2 max-w-[260px] truncate" title={term.term}>{term.term}</td>
                        <td className="px-4 py-2">{(term.adjP ?? term.pValue).toExponential(2)}</td>
                        <td className="px-4 py-2">{term.combinedScore?.toFixed(2) ?? '-'}</td>
                        <td className="px-4 py-2 max-w-[320px] truncate" title={(term.genes ?? []).join(', ')}>
                          {(term.genes ?? []).join(', ')}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </>
          )}

          {results.length === 0 && !loading && (
            <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md text-sm text-slate-500 dark:text-slate-400">
              Run enrichment to visualize pathway hits.
            </div>
          )}

          <div className="bg-white dark:bg-slate-900 border border-slate-200 dark:border-slate-800 rounded-xl p-6 shadow-md space-y-4">
            <div className="flex flex-col lg:flex-row lg:items-center lg:justify-between gap-4">
              <div>
                <h2 className="text-lg font-bold text-slate-900 dark:text-white flex items-center gap-2">
                  <Share2 className="w-5 h-5 text-cyan-600 dark:text-cyan-400" />
                  STRING PPI Network
                </h2>
                <p className="text-sm text-slate-500 dark:text-slate-400">
                  Builds a protein-protein interaction graph from high-expression genes.
                </p>
              </div>
              <Button onClick={runStringNetwork} isLoading={ppiLoading} className="justify-center">
                Generate PPI
              </Button>
            </div>

            {ppiStats && (
              <p className="text-xs text-slate-500 dark:text-slate-400">
                Filtered genes: {ppiStats.filtered} · Mapped to STRING: {ppiStats.mapped}
              </p>
            )}

            {ppiError && (
              <div className="p-3 rounded-lg bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-500/40 text-sm text-red-600 dark:text-red-300">
                {ppiError}
              </div>
            )}

            {ppiImage && (
              <div className="space-y-3">
                {ppiLink && (
                  <a
                    href={ppiLink}
                    target="_blank"
                    rel="noreferrer"
                    className="text-xs text-cyan-600 dark:text-cyan-400 underline"
                  >
                    Open interactive STRING view
                  </a>
                )}
                <div className="w-full overflow-x-auto rounded-lg border border-slate-200 dark:border-slate-800 bg-slate-50 dark:bg-slate-950/40 p-2">
                  <img src={ppiImage} alt="STRING network" className="w-full h-auto" />
                </div>
              </div>
            )}
          </div>
        </div>
      </div>
    </>
  );
};
