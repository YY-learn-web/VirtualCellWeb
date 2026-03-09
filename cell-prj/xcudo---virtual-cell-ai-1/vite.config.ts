import path from 'path';
import { defineConfig, loadEnv } from 'vite';
import react from '@vitejs/plugin-react';

export default defineConfig(({ mode }) => {
  const env = loadEnv(mode, '.', '');
  const apiProxyTarget = env.VITE_API_PROXY_TARGET || env.VITE_API_BASE_URL || 'http://127.0.0.1:8080';
  const proxy = {
    '/api': {
      target: apiProxyTarget,
      changeOrigin: true,
      secure: false,
    },
  };

  return {
    server: {
      port: 3000,
      host: '0.0.0.0',
      proxy,
    },
    preview: {
      port: 3000,
      host: '0.0.0.0',
      proxy,
    },
    plugins: [react()],
    define: {
      'process.env.API_KEY': JSON.stringify(env.GEMINI_API_KEY),
      'process.env.GEMINI_API_KEY': JSON.stringify(env.GEMINI_API_KEY),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '.'),
      },
    },
  };
});
