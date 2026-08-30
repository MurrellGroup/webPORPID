import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { Methods, type MethodsPage } from "./Methods";
import "./styles.css";

const page = (document.body.dataset.methodsPage ?? "index") as MethodsPage;
createRoot(document.getElementById("root")!).render(<StrictMode><Methods page={page} /></StrictMode>);
