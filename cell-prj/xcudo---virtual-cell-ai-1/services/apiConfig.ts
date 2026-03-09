const trimTrailingSlash = (value: string) => value.replace(/\/+$/, '');

const configuredApiBase = import.meta.env.VITE_API_BASE_URL?.trim();

const API_BASE_URL = configuredApiBase ? trimTrailingSlash(configuredApiBase) : '';

export const buildApiUrl = (path: string): string => {
  if (/^https?:\/\//i.test(path)) {
    return path;
  }

  const normalizedPath = path.startsWith('/') ? path : `/${path}`;
  return API_BASE_URL ? `${API_BASE_URL}${normalizedPath}` : normalizedPath;
};

export const mapNetworkError = (error: unknown): Error => {
  if (error instanceof TypeError) {
    return new Error(
      'Cannot reach the ASCEND backend. Start the FastAPI service on port 8080, or set VITE_API_BASE_URL to the backend origin.',
    );
  }

  if (error instanceof Error) {
    return error;
  }

  return new Error(String(error));
};

export { API_BASE_URL };
